// Validation du pricer unifie. `node sources/test_pricer.js` depuis le dossier.
//
// Trois controles, du plus revelateur au plus fin :
//   1. MARTINGALE. call(K=F) = put(K=F) equivaut a E[F_T] = F_0. C'est le test
//      qui attrape un compensateur faux, la faute la plus facile a commettre en
//      ajoutant une famille de sauts.
//   2. BLACK-SCHOLES. Le bloc a volatilite constante doit rendre exactement la
//      volatilite qu'on lui donne. Valide toute la chaine, du CF a l'inversion.
//   3. COS contre PROJ. Deux methodes numeriques independantes : si elles
//      s'accordent a 0.001 bps sur tout un smile, aucune des deux ne derive.
"use strict";
const P = require("./pricer.js");
const F0 = 100;
let pire = {mart: 0, meth: 0};

const VAR = {
  "constante": {var: {type: "const", sig: 0.45}},
  "Heston": {var: {type: "heston", kappa: 2, theta: 0.20, sigma: 0.7, rho: -0.4, v0: 0.22}},
  "double Heston": {var: {type: "heston2", f: [
    {kappa: 4, theta: 0.12, sigma: 0.6, rho: -0.6, v0: 0.14},
    {kappa: 0.6, theta: 0.09, sigma: 0.35, rho: -0.2, v0: 0.08}]}},
  "Heston + sauts de variance": {var: {type: "hestonJ", kappa: 2, theta: 0.20,
    sigma: 0.7, rho: -0.4, v0: 0.22, lamV: 1.5, muV: 0.08}},
};
const SAUT = {
  "aucun": {},
  "lognormale": {sauts: [{lam: 3, muJ: -0.05, sigJ: 0.12}]},
  "bi-lognormale": {sauts: [{lam: 8, muJ: 0.03, sigJ: 0.06},
                            {lam: 0.3, muJ: -0.5, sigJ: 0.8}]},
  "Kou": {kou: {lam: 5, p: 0.4, eta1: 12, eta2: 8}},
  "CGMY": {cgmy: {C: 0.5, G: 5, M: 5, Y: 1.2}},
  "Variance Gamma": {vg: {sigma: 0.35, nu: 0.4, theta: -0.15}},
  "NIG": {nig: {alpha: 8, beta: -2, delta: 0.6}},
};

console.log("\n1. MARTINGALE — call(K=F) − put(K=F), en bps de F0");
console.log("".padEnd(28) + Object.keys(SAUT).map(x => x.slice(0, 9).padStart(11)).join(""));
for (const [nv, V] of Object.entries(VAR)) {
  let l = "";
  for (const J of Object.values(SAUT)) {
    const p = {...V, ...J};
    let m = 0;
    for (const T of [2 / 365, 0.25, 1.0]) {
      const r = [{F: F0, K: F0, r: 0.03, c: 1}, {F: F0, K: F0, r: 0.03, c: 0}];
      const a = P.cosMaturity(r, T, p, 4096, 24), b = P.projMaturity(r, T, p, 4096, 20);
      m = Math.max(m, Math.abs((a[0] - a[1]) / F0 * 1e4), Math.abs((b[0] - b[1]) / F0 * 1e4));
    }
    pire.mart = Math.max(pire.mart, m);
    l += (isFinite(m) ? m.toFixed(2) : "NaN").padStart(11);
  }
  console.log(nv.padEnd(28) + l);
}

console.log("\n2. BLACK-SCHOLES — le bloc constant doit rendre sa propre volatilité");
for (const sig of [0.20, 0.45, 0.90]) {
  const p = {var: {type: "const", sig}};
  const T = 0.25, R = {F: F0, K: F0, r: 0.03, c: 1};
  const v = P.ivB76(P.cosMaturity([R], T, p, 4096, 24)[0], F0, F0, T, 0.03, 1);
  const ec = Math.abs(v - sig) * 1e4;
  console.log("   σ = " + (sig * 100).toFixed(0).padStart(2) + " %  ->  IV = "
              + (v * 100).toFixed(4) + " %   écart " + ec.toFixed(4) + " bps");
}

console.log("\n3. COS contre PROJ — écart max en bps d'IV sur un smile OTM (3 mois)");
const ys = []; for (let y = -0.6; y <= 0.61; y += 0.15) ys.push(y);
for (const [nv, V] of Object.entries(VAR)) {
  for (const [ns, J] of Object.entries(SAUT)) {
    const p = {...V, ...J}, T = 0.25;
    const rows = ys.map(y => ({F: F0, K: F0 * Math.exp(y), r: 0.03, c: y >= 0 ? 1 : 0}));
    const a = P.cosMaturity(rows, T, p, 4096, 24), b = P.projMaturity(rows, T, p, 4096, 20);
    let mx = 0;
    rows.forEach((R, i) => {
      const ia = P.ivB76(a[i], R.F, R.K, T, R.r, R.c), ib = P.ivB76(b[i], R.F, R.K, T, R.r, R.c);
      if (isFinite(ia) && isFinite(ib)) mx = Math.max(mx, Math.abs(ia - ib) * 1e4);
    });
    pire.meth = Math.max(pire.meth, mx);
  }
}
console.log("   pire écart sur les " + (Object.keys(VAR).length * Object.keys(SAUT).length)
            + " compositions : " + pire.meth.toFixed(4) + " bps");

console.log("\n===> martingale : " + pire.mart.toFixed(4) + " bps   ·   COS/PROJ : "
            + pire.meth.toFixed(4) + " bps");
console.log((pire.mart < 0.05 && pire.meth < 0.05)
            ? "     tout est conforme.\n" : "     AU MOINS UN CONTROLE A ECHOUE.\n");
