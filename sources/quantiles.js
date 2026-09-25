/* Quantiles de la loi diffusee, par inversion de la fonction caracteristique.

   POURQUOI PAS DE MONTE-CARLO. La methodologie demande 101 quantiles archives
   par date et par horizon (annexe C). Les tirer par simulation coute cher et,
   surtout, rend un rang qui SATURE : si le realise depasse toutes les
   trajectoires, le rang vaut exactement 1 et Anderson-Darling explose sur un
   log(0) -- on l'a vu en pratique. Ici les quantiles viennent de la meme
   fonction caracteristique que le pricer, donc :

     - toutes les compositions marchent sans echantillonneur dedie (CGMY et NIG
       demanderaient des representations en serie ou des subordinateurs) ;
     - aucun bruit de simulation : deux appels identiques rendent le meme rang ;
     - un rang ne sature jamais, la queue est evaluee analytiquement.

   METHODE. Densite par COS (Fang & Oosterlee) sur la fenetre de troncature
   donnee par les cumulants, puis cumul trapezoidal et inversion par
   interpolation lineaire. La fenetre est la meme que celle du pricer, ce qui
   evite d'avoir deux reglages a accorder.

   CE QU'ON DIFFUSE. X_T = log(F_T / F_0), le log-forward, qui est une
   martingale sous Q. Le facteur relatif d'un future est exp(X_T) et ne depend
   pas de l'echeance. Le COMPTANT, lui, porte la base : S_T/S_0 = exp(X_T + rT)
   avec r implicite Deribit -- c'est a l'appelant de l'ajouter. */
(function (root) {
  "use strict";
  const P = root.Pricer || root.BatesPricer || (typeof require === "function" ? require("./pricer.js") : null);

  // cumulants par differences finies sur logCF : suit le modele quel qu'il soit
  function cumulants(p, T) {
    const h = 1e-5;
    const l0 = P.logCF(0, T, p), lp = P.logCF(h, T, p), lm = P.logCF(-h, T, p);
    const c1 = (lp[1] - lm[1]) / (2 * h);
    const c2 = -(lp[0] - 2 * l0[0] + lm[0]) / (h * h);
    const c4 = P.c4Sauts ? P.c4Sauts(p, T) : 0;
    return {c1, c2, c4};
  }

  /* Densite de X_T sur une grille reguliere, par COS.
     Rend {x, f, a, b} ; f n'est pas renormalisee ici. */
  function densite(p, T, N, L) {
    N = N || 4096; L = L || 12;
    const {c1, c2, c4} = cumulants(p, T);
    if (!isFinite(c1) || !isFinite(c2) || c2 <= 0) return null;
    const demi = L * Math.sqrt(Math.abs(c2) + Math.sqrt(Math.abs(c4)));
    const a = c1 - demi, b = c1 + demi, ba = b - a;
    // coefficients A_k = (2/(b-a)) Re{ phi(k pi/(b-a)) exp(-i k pi a/(b-a)) }
    const A = new Float64Array(N);
    for (let k = 0; k < N; k++) {
      const u = k * Math.PI / ba;
      const l = P.logCF(u, T, p);
      const e = Math.exp(l[0]);
      const ang = l[1] - u * a;
      A[k] = (2 / ba) * e * Math.cos(ang);
    }
    A[0] *= 0.5;                                  // le prime de la somme
    const M = 2048;                               // points de la grille en x
    const x = new Float64Array(M), f = new Float64Array(M);
    for (let j = 0; j < M; j++) {
      const xx = a + (j + 0.5) * ba / M;
      x[j] = xx;
      let s = 0;
      const t = Math.PI * (xx - a) / ba;
      for (let k = 0; k < N; k++) s += A[k] * Math.cos(k * t);
      f[j] = s;
    }
    return {x, f, a, b};
  }

  /* Fonction de repartition sur la grille, monotone et normalisee a 1. */
  function repartition(p, T, N, L) {
    const d = densite(p, T, N, L);
    if (!d) return null;
    const M = d.x.length, dx = d.x[1] - d.x[0];
    const F = new Float64Array(M);
    let cum = 0;
    for (let j = 0; j < M; j++) {
      // la densite COS peut osciller legerement sous zero dans les queues :
      // on tronque, sinon F cesse d'etre monotone et l'inversion se casse
      const aire = Math.max(d.f[j], 0) * dx;
      // x[j] est le CENTRE de la cellule : F y vaut le cumul des cellules
      // precedentes plus la MOITIE de celle-ci. Attribuer le cumul complet au
      // centre decale toute la repartition d'une demi-maille -- un biais
      // systematique de dx/2, invisible sur la mediane mais qui domine dans
      // les queues, la ou la densite est plate.
      F[j] = cum + aire / 2;
      cum += aire;
    }
    // le total se lit maintenant sur le bord droit, pas au dernier centre
    const tot = cum;
    if (!(tot > 0.9 && tot < 1.1)) return null;   // fenetre inadaptee
    for (let j = 0; j < M; j++) F[j] /= tot;
    return {x: d.x, F};
  }

  /* Les 101 quantiles de la methodologie : 0.1 %, 1 %, 2 %, ..., 99 %, 99.9 %
     (annexe C). Rend les valeurs de X_T, pas du facteur. */
  function probas101() {
    const p = [0.001];
    for (let i = 1; i <= 99; i++) p.push(i / 100);
    p.push(0.999);
    return p;                                     // 101 valeurs
  }

  function quantiles(par, T, probs, N, L) {
    const R = repartition(par, T, N, L);
    if (!R) return null;
    probs = probs || probas101();
    const out = new Float64Array(probs.length);
    let j = 0;
    probs.forEach((q, i) => {
      while (j < R.F.length - 1 && R.F[j] < q) j++;
      if (j === 0) { out[i] = R.x[0]; return; }
      const F0 = R.F[j - 1], F1 = R.F[j];
      const w = F1 > F0 ? (q - F0) / (F1 - F0) : 0;
      out[i] = R.x[j - 1] + w * (R.x[j] - R.x[j - 1]);
    });
    return out;
  }

  /* Rang au sens de l'annexe C : Card(quantiles <= realise) / 101.
     `realise` est le log-rendement observe, X_T realise. */
  function rang(quants, realise) {
    if (!quants || !isFinite(realise)) return NaN;
    let n = 0;
    for (let i = 0; i < quants.length; i++) if (quants[i] <= realise) n++;
    return n / quants.length;
  }

  const api = {cumulants, densite, repartition, quantiles, probas101, rang};
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.Quantiles = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
