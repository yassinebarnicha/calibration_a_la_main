/* Tests statistiques de la methodologie de backtesting, facteur de risque.

   Reference : « Backtesting Methodology — Counterparty Credit Risk, EEPE /
   Internal Model Method », sections 3.2 a 3.4 et annexes B a G.

   TOUT PART DU RANG. Pour chaque date d'initialisation et chaque horizon, le
   modele archive 101 quantiles de sa distribution diffusee (0.1 %, 1 %, 2 %,
   ..., 99 %, 99.9 %) et le rang de la valeur realisee vaut

       rang = Card(quantiles <= valeur realisee) / 101          (annexe C)

   Si le modele decrit correctement le facteur, ces rangs sont i.i.d. uniformes
   sur [0, 1]. Tous les tests qui suivent mesurent un ecart a cette hypothese.

   DEUX TESTS TRAVAILLENT SUR LE RANG TRANSFORME. Kolmogorov-Smirnov et le chi2
   sur la variance appliquent d'abord z = Phi^-1(rang) et se comparent a la
   gaussienne standard -- c'est ce que prescrivent les annexes C et G, parce
   que l'implementation de reference compare a une gaussienne. Anderson-Darling
   et le test de signe travaillent sur le rang brut, uniforme.

   AUCUN GENERATEUR ALEATOIRE ICI. Les quantiles viennent de l'inversion de la
   fonction caracteristique, pas d'un Monte-Carlo : un rang ne peut donc jamais
   saturer a 0 ou a 1 faute de trajectoires. */
(function (root) {
  "use strict";

  // ---- fonctions speciales ------------------------------------------------
  // Lanczos g=7, n=9 : la meme que le pricer, pour ne pas avoir deux verites
  const LANCZOS = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                   771.32342877765313, -176.61502916214059, 12.507343278686905,
                   -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  function gamma(x) {
    if (x < 0.5) return Math.PI / (Math.sin(Math.PI * x) * gamma(1 - x));
    x -= 1;
    let a = LANCZOS[0];
    const t = x + 7.5;
    for (let i = 1; i < 9; i++) a += LANCZOS[i] / (x + i);
    return Math.sqrt(2 * Math.PI) * Math.pow(t, x + 0.5) * Math.exp(-t) * a;
  }
  function lgamma(x) {
    if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lgamma(1 - x);
    x -= 1;
    let a = LANCZOS[0];
    const t = x + 7.5;
    for (let i = 1; i < 9; i++) a += LANCZOS[i] / (x + i);
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  }

  // gamma incomplete reguliere P(a,x), serie puis fraction continue de Lentz
  function gammaP(a, x) {
    if (x <= 0) return 0;
    if (x < a + 1) {                      // serie : converge vite a gauche
      let ap = a, som = 1 / a, del = som;
      for (let n = 0; n < 500; n++) {
        ap++; del *= x / ap; som += del;
        if (Math.abs(del) < Math.abs(som) * 1e-15) break;
      }
      return som * Math.exp(-x + a * Math.log(x) - lgamma(a));
    }
    // fraction continue pour Q(a,x), puis P = 1 - Q
    let b = x + 1 - a, c = 1e300, d = 1 / b, h = d;
    for (let i = 1; i < 500; i++) {
      const an = -i * (i - a);
      b += 2; d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d;
      const del = d * c; h *= del;
      if (Math.abs(del - 1) < 1e-15) break;
    }
    return 1 - Math.exp(-x + a * Math.log(x) - lgamma(a)) * h;
  }
  const chi2cdf = (x, df) => (df > 0 && x > 0) ? gammaP(df / 2, x / 2) : 0;

  // erf par la meme gamma incomplete : une seule primitive a faire confiance
  const erf = x => (x >= 0 ? 1 : -1) * gammaP(0.5, x * x);
  const normcdf = x => 0.5 * (1 + erf(x / Math.SQRT2));
  // Acklam, raffine par une iteration de Halley : ~1e-15
  function norminv(p) {
    if (!(p > 0 && p < 1)) return p <= 0 ? -Infinity : Infinity;
    const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
    const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
               6.680131188771972e+01, -1.328068155288572e+01];
    const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
               -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
    const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
               3.754408661907416e+00];
    const pl = 0.02425;
    let q, r, x;
    if (p < pl) {
      q = Math.sqrt(-2 * Math.log(p));
      x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    } else if (p <= 1 - pl) {
      q = p - 0.5; r = q * q;
      x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
          (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    } else {
      q = Math.sqrt(-2 * Math.log(1 - p));
      x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
    const e = normcdf(x) - p;
    const u = e * Math.sqrt(2 * Math.PI) * Math.exp(x * x / 2);
    return x - u / (1 + x * u / 2);
  }

  // binomiale cumulee P(X <= k), par la beta incomplete pour rester stable
  function binocdf(k, n, p) {
    if (k < 0) return 0;
    if (k >= n) return 1;
    let som = 0;
    // n reste petit ici (nombre d'observations) : somme directe en log
    for (let i = 0; i <= k; i++) {
      som += Math.exp(lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1)
                      + i * Math.log(p) + (n - i) * Math.log(1 - p));
    }
    return Math.min(1, som);
  }

  // ---- zones RAG ----------------------------------------------------------
  // p-valeur : rouge < 0.01 %, ambre < 5 %, vert au-dela (annexes C a G)
  const zoneP = p => (p < 0.0001 ? "rouge" : p < 0.05 ? "ambre" : "vert");
  // P_c des feux : vert dans [2.5 %, 97.5 %], rouge hors [0.005 %, 99.995 %]
  const zonePc = pc => (pc < 0.00005 || pc > 0.99995) ? "rouge"
                     : (pc < 0.025 || pc > 0.975) ? "ambre" : "vert";

  // ---- 1. feux tricolores (section 3.2, annexe B) -------------------------
  // Intervalles de quantiles NON cumulatifs, par defaut des deciles. Pour
  // chaque intervalle : combien de rangs y tombent, contre combien on en
  // attendait, et la probabilite binomiale cumulee d'un tel comptage.
  function feux(rangs, nb) {
    nb = nb || 10;
    const T = rangs.length, p = 1 / nb, attendu = T * p;
    const cases = [];
    for (let i = 0; i < nb; i++) {
      const lo = i / nb, hi = (i + 1) / nb;
      // dernier intervalle ferme a droite, pour que 1.0 soit compte
      const n = rangs.filter(r => r >= lo && (i === nb - 1 ? r <= hi : r < hi)).length;
      // Convention : P(X <= n), binomiale cumulee usuelle. La source ecrit une
      // somme demarrant a k = 1, que sa propre note de reconstruction signale
      // comme douteuse ; on retient la forme standard.
      const pc = binocdf(n, T, p);
      cases.push({lo, hi, attendu, realise: n,
                  ratio: attendu > 0 ? n / attendu : NaN, pc, zone: zonePc(pc)});
    }
    const pire = cases.reduce((a, c) =>
      ({rouge: 3, ambre: 2, vert: 1}[c.zone] > {rouge: 3, ambre: 2, vert: 1}[a] ? c.zone : a), "vert");
    return {cases, zone: pire, T, nb};
  }

  // ---- 2. Kolmogorov-Smirnov (section 3.3.1, annexe C) --------------------
  // Les rangs sont transformes par Phi^-1 et compares a la gaussienne standard.
  function ks(rangs) {
    const n = rangs.length;
    if (n < 3) return {n, stat: NaN, p: NaN, zone: "—"};
    const z = rangs.map(r => norminv(Math.min(Math.max(r, 1 / (2 * 101)), 1 - 1 / (2 * 101))))
                   .sort((a, b) => a - b);
    let D = 0;
    for (let i = 0; i < n; i++) {
      const F = normcdf(z[i]);
      D = Math.max(D, (i + 1) / n - F, F - i / n);
    }
    // Kolmogorov asymptotique, correction de Stephens pour n fini
    const lam = (Math.sqrt(n) + 0.12 + 0.11 / Math.sqrt(n)) * D;
    let som = 0;
    for (let k = 1; k <= 100; k++) som += Math.pow(-1, k - 1) * Math.exp(-2 * k * k * lam * lam);
    const p = Math.min(1, Math.max(0, 2 * som));
    return {n, stat: D, p, zone: zoneP(p)};
  }

  // ---- 3. Anderson-Darling (section 3.3.2, annexe D) ----------------------
  // Sur les rangs bruts, contre l'uniforme. p-valeur par Marsaglia.
  function adinf(z) {
    if (z < 2) {
      return Math.exp(-1.2337141 / z) / Math.sqrt(z) *
        (2.00012 + (0.247105 - (0.0649821 - (0.0347962 - (0.011672 - 0.00168691 * z)
          * z) * z) * z) * z);
    }
    return Math.exp(-Math.exp(1.0776 - (2.30695 - (0.43424 - (0.082433 -
      (0.008056 - 0.0003146 * z) * z) * z) * z) * z));
  }
  function errfix(n, x) {
    const g1 = x => Math.sqrt(x) * (1 - x) * (49 * x - 102);
    const g2 = x => -0.00022633 + (6.54034 - (14.6538 - (14.458 - (8.259 - 1.91864 * x)
      * x) * x) * x) * x;
    const g3 = x => -130.2137 + (745.2337 - (1705.091 - (1950.646 - (1116.360 -
      255.7844 * x) * x) * x) * x) * x;
    const c = 0.01265 + 0.1757 / n;
    if (x < c) { const t = x / c; return (0.0037 / (n * n * n) + 0.00078 / (n * n) + 0.00006 / n) * g1(t); }
    if (x < 0.8) { const t = (x - c) / (0.8 - c); return (0.04213 / n + 0.01365 / (n * n)) * g2(t); }
    return g3(x) / n;
  }
  function ad(rangs) {
    const n = rangs.length;
    if (n < 3) return {n, stat: NaN, p: NaN, zone: "—"};
    const x = rangs.map(r => Math.min(Math.max(r, 1e-12), 1 - 1e-12)).sort((a, b) => a - b);
    let som = 0;
    for (let k = 1; k <= n; k++) som += (2 * k - 1) * Math.log(x[k - 1] * (1 - x[n - k]));
    const A2 = -n - som / n;
    const p = Math.min(1, Math.max(0, 1 - adinf(A2) - errfix(n, adinf(A2))));
    return {n, stat: A2, p, zone: zoneP(p)};
  }

  // ---- 4. Ljung-Box (section 3.4.1, annexe E) -----------------------------
  // m se deduit du recouvrement : horizon / pas d'initialisation, moins un.
  function ljungBox(rangs, m) {
    const N = rangs.length;
    m = Math.max(1, Math.min(m || 1, Math.floor(N / 4)));
    if (N < 8) return {N, m, stat: NaN, p: NaN, zone: "—"};
    const mu = rangs.reduce((a, b) => a + b, 0) / N;
    const c = rangs.map(r => r - mu);
    const c0 = c.reduce((a, b) => a + b * b, 0);
    let Q = 0;
    const acf = [];
    for (let h = 1; h <= m; h++) {
      let s = 0;
      for (let t = 0; t + h < N; t++) s += c[t] * c[t + h];
      const rho = c0 > 0 ? s / c0 : 0;
      acf.push(rho);
      Q += rho * rho / (N - h);
    }
    Q *= N * (N + 2);
    const p = 1 - chi2cdf(Q, m);
    return {N, m, stat: Q, p, zone: zoneP(p), acf};
  }

  // ---- 5. test de signe bilateral sur la mediane (3.4.2, annexe F) --------
  // L'echantillon est centre en retranchant 0.5 ; on compte au-dessus et en
  // dessous. Exact par la binomiale, gaussien au-dela de 100 points.
  function signe(rangs) {
    const N = rangs.length;
    if (N < 5) return {N, p: NaN, zone: "—"};
    const pos = rangs.filter(r => r > 0.5).length;
    const neg = rangs.filter(r => r < 0.5).length;
    let p, stat;
    if (N < 100) {
      stat = Math.min(pos, neg);
      p = Math.min(1, 2 * binocdf(stat, pos + neg, 0.5));
    } else {
      const sg = (pos - neg - Math.sign(pos - neg)) / Math.sqrt(N);
      stat = sg;
      p = 2 * normcdf(-Math.abs(sg));
    }
    return {N, pos, neg, stat, p: Math.min(1, p), zone: zoneP(p),
            // au-dessus de 0.5 : le modele SOUS-estime le facteur (3.4.2)
            sens: pos > neg ? "sous-estimation" : pos < neg ? "surestimation" : "neutre"};
  }

  // ---- 6. chi2 bilateral sur la variance (3.4.3, annexe G) ---------------
  // Sur les rangs transformes par Phi^-1 : la variance de reference vaut 1,
  // ce qui leve au passage l'ambiguite s/v de la source.
  function variance(rangs) {
    const N = rangs.length;
    if (N < 5) return {N, p: NaN, zone: "—"};
    const z = rangs.map(r => norminv(Math.min(Math.max(r, 1 / (2 * 101)), 1 - 1 / (2 * 101))));
    const mu = z.reduce((a, b) => a + b, 0) / N;
    const s2 = z.reduce((a, b) => a + (b - mu) * (b - mu), 0) / (N - 1);
    const T = (N - 1) * s2;                       // v = 1 apres transformation
    const c = chi2cdf(T, N - 1);
    const p = Math.min(1, 2 * Math.min(c, 1 - c));
    // le ratio rapporte dans le tableau reste sur les rangs bruts : 1/12
    const mub = rangs.reduce((a, b) => a + b, 0) / N;
    const vb = rangs.reduce((a, b) => a + (b - mub) * (b - mub), 0) / (N - 1);
    return {N, stat: T, p, zone: zoneP(p), varObs: vb, ratio: vb / (1 / 12),
            sens: vb > 1 / 12 ? "volatilité sous-estimée" : "volatilité surestimée"};
  }

  // ---- synthese pour un couple (facteur de risque, horizon) --------------
  function analyser(rangs, opts) {
    opts = opts || {};
    const r = rangs.filter(x => isFinite(x));
    const moy = r.length ? r.reduce((a, b) => a + b, 0) / r.length : NaN;
    return {
      n: r.length,
      rangs: r,
      moyenne: moy,
      moyenneRatio: moy / 0.5,                    // « observé / attendu » du tableau
      feux: feux(r, opts.intervalles || 10),
      ks: ks(r),
      ad: ad(r),
      ljungBox: opts.recouvrement ? ljungBox(r, opts.lags) : null,
      signe: signe(r),
      variance: variance(r),
    };
  }

  const api = {gamma, lgamma, chi2cdf, normcdf, norminv, binocdf,
               feux, ks, ad, ljungBox, signe, variance, analyser,
               zoneP, zonePc};
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.Backtest = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
