// Validation des quantiles par inversion de Fourier.
//
// Le controle decisif est Black-Scholes : avec le bloc a volatilite constante
// et aucun saut, X_T = log(F_T/F_0) suit exactement une gaussienne de moyenne
// -sigma^2 T/2 et d'ecart-type sigma sqrt(T). On a donc une reference ANALYTIQUE
// pour les 101 quantiles, pas une comparaison a un autre calcul approche.
"use strict";
const Q = require("./quantiles.js");
const B = require("./backtest.js");
let ko = 0;
const ok = (c, m) => { if (!c) ko++; console.log((c ? "  ok   " : "  ECHEC") + " " + m); };

console.log("\n1. BLACK-SCHOLES : quantiles contre la gaussienne exacte");
for (const [sig, T] of [[0.45, 1 / 365], [0.45, 30 / 365], [0.80, 0.5], [0.25, 1.0]]) {
  const p = {var: {type: "const", sig}};
  const probs = Q.probas101();
  const num = Q.quantiles(p, T, probs);
  if (!num) { ok(false, "sigma=" + sig + " T=" + T.toFixed(4) + " : pas de quantiles"); continue; }
  const mu = -0.5 * sig * sig * T, sd = sig * Math.sqrt(T);
  let pire = 0, pireP = 0;
  probs.forEach((q, i) => {
    const exact = mu + sd * B.norminv(q);
    const e = Math.abs(num[i] - exact) / sd;         // ecart en ecarts-types
    if (e > pire) { pire = e; pireP = q; }
  });
  ok(pire < 2e-3, "sigma=" + (sig * 100).toFixed(0) + "% T=" + (T * 365).toFixed(0)
     + "j : ecart max " + (pire * 1e4).toFixed(2) + " pb d'ecart-type (au quantile "
     + (pireP * 100).toFixed(1) + "%)");
}

console.log("\n2. LA REPARTITION EST UNE VRAIE LOI");
const jeux = {
  "Heston + bi-lognormale": {var: {type: "heston", kappa: 3, theta: .16, sigma: .8, rho: -.3, v0: .16},
    sauts: [{lam: 10, muJ: .03, sigJ: .06}, {lam: .2, muJ: -.4, sigJ: .7}]},
  "double Heston + CGMY": {var: {type: "heston2", f: [
      {kappa: 6, theta: .12, sigma: .7, rho: -.6, v0: .14},
      {kappa: .6, theta: .12, sigma: .35, rho: -.2, v0: .08}]},
    cgmy: {C: .3, G: 6, M: 6, Y: 1.25}},
  "Heston + sauts de variance + Kou": {var: {type: "hestonJ", kappa: 3, theta: .16,
      sigma: .8, rho: -.3, v0: .16, lamV: 1.5, muV: .08},
    kou: {lam: 8, p: .4, eta1: 12, eta2: 8}},
  "constante + NIG": {var: {type: "const", sig: .5}, nig: {alpha: 10, beta: -3, delta: .5}},
};
for (const [nom, p] of Object.entries(jeux)) {
  const T = 21 / 365;
  const R = Q.repartition(p, T);
  const qs = Q.quantiles(p, T);
  if (!R || !qs) { ok(false, nom + " : echec"); continue; }
  let monotone = true;
  for (let i = 1; i < qs.length; i++) if (qs[i] < qs[i - 1]) monotone = false;
  // E[e^X] doit valoir 1 : c'est la martingale, relue par la densite
  const dx = R.x[1] - R.x[0];
  let esp = 0;
  for (let j = 1; j < R.x.length; j++) esp += Math.exp(R.x[j]) * (R.F[j] - R.F[j - 1]);
  const med = qs[50];                               // le quantile 50 %
  ok(monotone && Math.abs(esp - 1) < 5e-3,
     nom.padEnd(34) + " monotone, E[e^X] = " + esp.toFixed(5)
     + ", mediane = " + (med * 100).toFixed(2) + " %");
}

console.log("\n3. LE RANG SE COMPORTE COMME IL DOIT");
const p3 = jeux["Heston + bi-lognormale"], T3 = 21 / 365;
const q3 = Q.quantiles(p3, T3);
ok(q3.length === 101, "101 quantiles archives, comme l'annexe C");
ok(Math.abs(Q.rang(q3, q3[50]) - 0.5) < 0.02,
   "rang de la mediane = " + Q.rang(q3, q3[50]).toFixed(3));
ok(Q.rang(q3, q3[0] - 1) === 0, "tres en dessous -> rang 0");
ok(Q.rang(q3, q3[100] + 1) === 1, "tres au-dessus -> rang 1");
// et surtout : un realise plausible ne sature pas
ok(Q.rang(q3, q3[98]) > 0.9 && Q.rang(q3, q3[98]) < 1,
   "au quantile 99 % : rang = " + Q.rang(q3, q3[98]).toFixed(3) + " (ne sature pas)");

console.log("\n4. REPRODUCTIBILITE : deux appels, le meme rang");
const a1 = Q.rang(Q.quantiles(p3, T3), 0.03);
const a2 = Q.rang(Q.quantiles(p3, T3), 0.03);
ok(a1 === a2, "rang identique a l'identique : " + a1.toFixed(4)
   + "  (un Monte-Carlo en aurait donne deux)");

console.log("\n===> " + (ko === 0 ? "tout est conforme." : ko + " CONTROLE(S) EN ECHEC.") + "\n");
process.exit(ko ? 1 : 0);
