// Validation du module statistique. `node sources/test_backtest.js`
//
// Le controle central est la TABLE DE L'ANNEXE B : la methodologie y publie
// 136 observations, leur repartition en dix deciles et les P_c associes. Si
// nos P_c reproduisent les siens, la formule est bonne ET la convention de
// sommation est tranchee -- la source ecrit une somme demarrant a k = 1, que
// sa propre note de reconstruction signale comme douteuse.
"use strict";
const B = require("./backtest.js");
let ko = 0;
const ok = (c, m) => { if (!c) ko++; console.log((c ? "  ok   " : "  ECHEC") + " " + m); };
const pres = (a, b, tol) => Math.abs(a - b) <= tol;

console.log("\n1. FONCTIONS SPECIALES");
ok(pres(B.chi2cdf(3.8415, 1), 0.95, 1e-4), "chi2cdf(3.8415, 1) = 0.95 -> " + B.chi2cdf(3.8415, 1).toFixed(6));
ok(pres(B.chi2cdf(5.9915, 2), 0.95, 1e-4), "chi2cdf(5.9915, 2) = 0.95 -> " + B.chi2cdf(5.9915, 2).toFixed(6));
ok(pres(B.chi2cdf(18.307, 10), 0.95, 1e-4), "chi2cdf(18.307, 10) = 0.95 -> " + B.chi2cdf(18.307, 10).toFixed(6));
ok(pres(B.normcdf(1.959964), 0.975, 1e-7), "normcdf(1.959964) = 0.975 -> " + B.normcdf(1.959964).toFixed(9));
ok(pres(B.norminv(0.975), 1.959964, 1e-6), "norminv(0.975) = 1.959964 -> " + B.norminv(0.975).toFixed(9));
ok(pres(B.norminv(0.001), -3.090232, 1e-6), "norminv(0.001) = -3.090232 -> " + B.norminv(0.001).toFixed(9));
ok(pres(B.binocdf(5, 10, 0.5), 0.623046875, 1e-12), "binocdf(5,10,0.5) = 0.623046875 -> " + B.binocdf(5, 10, 0.5));

console.log("\n2. TABLE DE L'ANNEXE B  (T = 136, p = 0.1, attendu 13.6 par decile)");
const ANNEXE_B = [
  [12, 0.88, 39.01], [13, 0.96, 50.40], [18, 1.32, 91.49], [9, 0.66, 11.68],
  [13, 0.96, 50.40], [11, 0.81, 28.26], [18, 1.32, 91.49], [18, 1.32, 91.49],
  [16, 1.18, 80.00], [8, 0.59, 6.51]];
const T = ANNEXE_B.reduce((a, x) => a + x[0], 0);
ok(T === 136, "les comptages somment a 136 -> " + T);
console.log("  décile   réalisé   ratio publié / calculé    P_c publié / calculé");
let pireEcart = 0;
ANNEXE_B.forEach(([n, ratio, pc], i) => {
  const r = n / (136 * 0.1);
  const p = B.binocdf(n, 136, 0.1) * 100;
  pireEcart = Math.max(pireEcart, Math.abs(p - pc));
  console.log("  " + String(i * 10).padStart(3) + "-" + String((i + 1) * 10).padEnd(4)
    + String(n).padStart(6) + "   " + ratio.toFixed(2) + " / " + r.toFixed(2)
    + "        " + pc.toFixed(2).padStart(6) + " % / " + p.toFixed(2).padStart(6) + " %");
});
ok(pireEcart < 0.02, "P_c reproduits a " + pireEcart.toFixed(4) + " point de pourcentage pres");
ok(B.zonePc(0.3901) === "vert" && B.zonePc(0.0651) === "vert",
   "les dix P_c publies tombent bien en zone verte");

console.log("\n3. UNIFORMITE : un echantillon uniforme ne doit pas etre rejete");
const uni = [];
for (let i = 1; i <= 200; i++) uni.push((i - 0.5) / 200);   // grille reguliere
ok(B.ks(uni).p > 0.05, "KS sur grille uniforme : p = " + B.ks(uni).p.toFixed(4));
ok(B.ad(uni).p > 0.05, "AD sur grille uniforme : p = " + B.ad(uni).p.toFixed(4));
ok(B.signe(uni).p > 0.05, "signe : p = " + B.signe(uni).p.toFixed(4));
ok(B.variance(uni).p > 0.05, "variance : p = " + B.variance(uni).p.toFixed(4)
   + ", ratio = " + B.variance(uni).ratio.toFixed(3));

console.log("\n4. NON-UNIFORMITE : un echantillon biaise doit etre rejete");
// tout dans la moitie haute : le modele sous-estime le facteur
const haut = uni.map(u => 0.5 + u / 2);
ok(B.ks(haut).p < 1e-4, "KS sur echantillon decale : p = " + B.ks(haut).p.toExponential(2));
ok(B.ad(haut).p < 1e-4, "AD : p = " + B.ad(haut).p.toExponential(2));
ok(B.signe(haut).p < 1e-4, "signe : p = " + B.signe(haut).p.toExponential(2)
   + " (" + B.signe(haut).sens + ")");
// tout au centre : le modele surestime la dispersion
const centre = uni.map(u => 0.35 + u * 0.3);
ok(B.variance(centre).p < 1e-4, "variance sur echantillon resserre : p = "
   + B.variance(centre).p.toExponential(2) + " (" + B.variance(centre).sens + ")");
ok(B.feux(haut).zone === "rouge", "feux sur echantillon decale : zone " + B.feux(haut).zone);

console.log("\n5. LJUNG-BOX : independance contre autocorrelation");
let x = 0.5;
const iid = [], ar1 = [];
let g = 12345;
const rnd = () => { g = (g * 1103515245 + 12345) % 2147483648; return g / 2147483648; };
for (let i = 0; i < 200; i++) {
  iid.push(rnd());
  x = 0.85 * (x - 0.5) + 0.15 * (rnd() - 0.5) + 0.5;     // fortement correle
  ar1.push(Math.min(Math.max(x, 0.001), 0.999));
}
ok(B.ljungBox(iid, 5).p > 0.05, "iid : p = " + B.ljungBox(iid, 5).p.toFixed(4));
ok(B.ljungBox(ar1, 5).p < 1e-4, "AR(1) : p = " + B.ljungBox(ar1, 5).p.toExponential(2)
   + ", rho1 = " + B.ljungBox(ar1, 5).acf[0].toFixed(3));

console.log("\n6. SYNTHESE : analyser() assemble le tout");
const a = B.analyser(uni, {recouvrement: true, lags: 3});
ok(a.n === 200 && a.feux.cases.length === 10 && a.ljungBox !== null,
   "n = " + a.n + ", " + a.feux.cases.length + " deciles, moyenne/attendu = "
   + (a.moyenneRatio * 100).toFixed(1) + " %");

console.log("\n===> " + (ko === 0 ? "tout est conforme." : ko + " CONTROLE(S) EN ECHEC.") + "\n");
process.exit(ko ? 1 : 0);
