/* Onglet Backtest — facteur de risque, selon la methodologie Natixis.

   CE QU'ON BACKTESTE. Uniquement des FACTEURS DE RISQUE, pas de portefeuille :
   les futures cotes, contrat par contrat, et le comptant.

   LE COMPTANT EST DEDUIT DE r. Sous Q, c'est le FUTURE qui est une martingale,
   pas le spot : S_t e^{-rt} l'est. Le rendement comptant comparable a la loi
   diffusee est donc log(S_{t+h}/S_t) - r*h, ou r est le taux implicite Deribit
   lu dans l'instantane, r = log(F/X)/T, moyenne au-dela de 20 jours -- en
   deca, T minuscule amplifie par 365 le moindre decalage d'horodatage et rend
   un taux sans aucun sens economique.

   LE RANG suit l'annexe C : 101 quantiles archives, rang = Card(q <= realise)
   / 101. Les quantiles viennent de l'inversion de Fourier (quantiles.js), pas
   d'un Monte-Carlo, donc le rang ne sature jamais et deux passages rendent le
   meme resultat.

   DEUX PLAGES INDEPENDANTES, comme demande :
     - la plage des PARAMETRES : de quels jours on tire les jeux calibres ;
     - la plage du BACKTESTING : sur quelle fenetre on observe.
   Elles ne se recouvrent pas forcement : on peut calibrer sur aout et
   observer septembre. */
(function (root) {
  "use strict";
  const Q = root.Quantiles, B = root.Backtest;

  // ---------- utilitaires de dates -----------------------------------------
  const jms = s => Date.parse(s.slice(0, 10) + "T00:00:00Z");
  const jplus = (s, n) => new Date(jms(s) + n * 864e5).toISOString().slice(0, 10);
  const ecartJours = (a, b) => Math.round((jms(b) - jms(a)) / 864e5);

  /* Taux implicite d'un instantane : moyenne au-dela de 20 jours. */
  function tauxImplicite(d) {
    const X = +d.index, vus = new Set(), rs = [];
    d.rows.forEach(r => {
      if (vus.has(r.m)) return;
      vus.add(r.m);
      if (r.T * 365 >= 20) rs.push(Math.log(r.F / X) / r.T);
    });
    return rs.length ? rs.reduce((a, b) => a + b, 0) / rs.length : 0;
  }

  /* Prix a terme par echeance dans un instantane deja converti en lignes. */
  function forwards(d) {
    const out = {}, vus = new Set();
    d.rows.forEach(r => { if (!vus.has(r.m)) { vus.add(r.m); out[d.mats[r.m].date] = r.F; } });
    return out;
  }

  // ---------- le moteur -----------------------------------------------------
  /* opts :
       ancres      [dates]      jours d'ou partent les previsions
       params      (jour) -> {} jeu de parametres a utiliser ce jour-la
       horizons    [entiers]    en jours calendaires
       facteurs    [{type:"future", ech} | {type:"spot"}]
       tolerance   entier       ecart accepte sur la date cible
       charger     (jour) -> D  acces aux instantanes (peut etre asynchrone)
     Rend {cles: [...], par: {cle: {horizon: [rangs]}}, detail: [...]} */
  async function lancer(opts, avance) {
    const {ancres, params, horizons, facteurs, charger} = opts;
    const tol = opts.tolerance == null ? 1 : opts.tolerance;
    const dispo = new Set(opts.disponibles || []);
    const par = {}, detail = [];
    const cle = f => f.type === "spot" ? "comptant" : "future " + f.ech;
    facteurs.forEach(f => { par[cle(f)] = {}; horizons.forEach(h => { par[cle(f)][h] = []; }); });

    // les quantiles ne dependent que (parametres, horizon) : on les partage
    const cacheQ = new Map();
    let fait = 0;
    for (const a of ancres) {
      const p = params(a);
      if (!p) continue;
      let dA; try { dA = await charger(a); } catch (e) { continue; }
      if (!dA) continue;
      const FA = forwards(dA), rA = tauxImplicite(dA), XA = +dA.index;

      for (const h of horizons) {
        // la date cible, a +/- tolerance : nos instantanes ont des trous
        let cible = null;
        for (let d = 0; d <= tol && !cible; d++) {
          for (const sg of (d === 0 ? [0] : [1, -1])) {
            const c = jplus(a, h + sg * d);
            if (dispo.has(c) && c > a) { cible = c; break; }
          }
        }
        if (!cible) continue;
        const hr = ecartJours(a, cible);            // l'horizon REEL
        let dB; try { dB = await charger(cible); } catch (e) { continue; }
        if (!dB) continue;
        const FB = forwards(dB), XB = +dB.index;

        const ck = JSON.stringify(p) + "|" + hr;
        let quants = cacheQ.get(ck);
        if (quants === undefined) {
          quants = Q.quantiles(p, hr / 365) || null;
          cacheQ.set(ck, quants);
        }
        if (!quants) continue;

        for (const f of facteurs) {
          let x = NaN, base = NaN, obs = NaN;
          if (f.type === "future") {
            if (!(f.ech in FA) || !(f.ech in FB)) continue;
            base = FA[f.ech]; obs = FB[f.ech];
            x = Math.log(obs / base);               // martingale : rien a retirer
          } else {
            base = XA; obs = XB;
            // le comptant porte la base : on la retire pour revenir a la martingale
            x = Math.log(obs / base) - rA * hr / 365;
          }
          if (!isFinite(x)) continue;
          const rg = Q.rang(quants, x);
          par[cle(f)][h].push(rg);
          detail.push({ancre: a, cible, horizon: h, horizonReel: hr, facteur: cle(f),
                       base, obs, logRendement: x, rang: rg});
        }
      }
      if (avance) avance(++fait, ancres.length);
      // on rend la main au navigateur : sinon la page gele sur une longue plage
      await new Promise(r => setTimeout(r, 0));
    }
    return {cles: facteurs.map(cle), par, detail};
  }

  /* m de Ljung-Box : avec des ancres espacees de `pas` jours et un horizon h,
     ceil(h/pas) - 1 observations consecutives partagent du temps (annexe E). */
  const lagsPour = (h, pas) => Math.max(1, Math.ceil(h / Math.max(pas, 1)) - 1);

  const api = {lancer, tauxImplicite, forwards, lagsPour, jplus, ecartJours};
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.MoteurBacktest = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
