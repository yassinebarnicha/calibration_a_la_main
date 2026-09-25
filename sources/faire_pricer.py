# -*- coding: utf-8 -*-
u"""Fabrique le pricer unifie de calibration_a_la_main.

On part du pricer CGMY -- qui porte deja Heston, les composantes lognormales et
CGMY -- et on lui ajoute Kou, Variance Gamma et NIG. COS et PROJ ne sont pas
touches : toutes les familles passent par le MEME logCF, donc par la meme
machinerie de prix. C'est tout l'interet de rester dans ce cadre.

Chaque famille apporte deux choses et deux seulement :
  psi(u)  l'exposant caracteristique de sa partie sauts,
  omega   son compensateur, psi(-i), reel, qui garantit E[F_T] = F_0.
et une formule fermee de quatrieme cumulant, qui ne sert qu'a dimensionner la
grille -- une mauvaise valeur n'y fait pas diverger le prix, elle degrade
silencieusement sa precision.
"""
import io
import shutil
from pathlib import Path

RAC = Path(r"C:/Users/busta/stage/Bates_calibrator")
SRC = RAC / "cgmy_a_la_main" / "sources" / "pricer.js"
DST = RAC / "calibration_a_la_main" / "sources"
DST.mkdir(parents=True, exist_ok=True)

s = io.open(SRC, encoding="utf-8").read()

ENTETE = u"""/* Pricer unifie : Heston + une famille de sauts au choix, COS et PROJ.

   Toutes les familles passent par le meme logCF, donc par la meme machinerie de
   prix -- c'est ce qui permet de jongler entre modeles sans rien changer
   d'autre. Une famille apporte son exposant psi(u) et son compensateur
   omega = psi(-i), reel, qui impose E[F_T] = F_0.

     p.sauts = [{lam, muJ, sigJ}, ...]   lognormales : Merton (1), bi-log (2)
     p.cgmy  = {C, G, M, Y}              CGMY, activite eventuellement infinie
     p.kou   = {lam, p, eta1, eta2}      double exponentielle asymetrique
     p.vg    = {sigma, nu, theta}        Variance Gamma
     p.nig   = {alpha, beta, delta}      Normal Inverse Gaussian

   Sans aucune de ces clefs : Heston pur. Avec sigma = 0 dans le bloc Heston :
   Levy pur. */
"""
s = s[s.index("(function (root)"):]
s = ENTETE + s

# ---------------------------------------------------------------- 1. les psi
NOUVEAU = u'''
  // ---- Kou : double exponentielle asymetrique ---------------------------
  // f(y) = p eta1 e^{-eta1 y} 1_{y>0} + (1-p) eta2 e^{eta2 y} 1_{y<0}
  // phi(u) = p eta1/(eta1 - iu) + (1-p) eta2/(eta2 + iu)
  // eta1 > 1 est OBLIGATOIRE : sinon E[e^Y] diverge et le compensateur n'existe
  // pas -- c'est la contrainte que les bornes doivent porter.
  function psiKou(u, j) {
    const a = cdiv([j.p * j.eta1, 0], [j.eta1, -u]);
    const b = cdiv([(1 - j.p) * j.eta2, 0], [j.eta2, u]);
    return [j.lam * (a[0] + b[0] - 1), j.lam * (a[1] + b[1])];
  }
  function compensateurKou(j) {
    return j.lam * (j.p * j.eta1 / (j.eta1 - 1) + (1 - j.p) * j.eta2 / (j.eta2 + 1) - 1);
  }

  // ---- Variance Gamma ---------------------------------------------------
  // psi(u) = -(1/nu) log(1 - i theta nu u + sigma^2 nu u^2 / 2)
  // Activite infinie, variation finie. 1 - theta nu - sigma^2 nu / 2 > 0 est
  // la condition d'existence du compensateur.
  function psiVG(u, j) {
    const z = [1 + 0.5 * j.sigma * j.sigma * j.nu * u * u, -j.theta * j.nu * u];
    const l = clog(z);
    return [-l[0] / j.nu, -l[1] / j.nu];
  }
  function compensateurVG(j) {
    const a = 1 - j.theta * j.nu - 0.5 * j.sigma * j.sigma * j.nu;
    return a > 0 ? -Math.log(a) / j.nu : NaN;
  }

  // ---- NIG --------------------------------------------------------------
  // psi(u) = -delta ( sqrt(alpha^2 - (beta + iu)^2) - sqrt(alpha^2 - beta^2) )
  // alpha > |beta + 1| pour que le compensateur existe.
  function psiNIG(u, j) {
    const r0 = Math.sqrt(Math.max(j.alpha * j.alpha - j.beta * j.beta, 0));
    // alpha^2 - (beta + iu)^2 = alpha^2 - beta^2 + u^2 - 2 i beta u
    const z = csqrt([j.alpha * j.alpha - j.beta * j.beta + u * u, -2 * j.beta * u]);
    return [-j.delta * (z[0] - r0), -j.delta * z[1]];
  }
  function compensateurNIG(j) {
    const a = j.alpha * j.alpha - (j.beta + 1) * (j.beta + 1);
    const b = j.alpha * j.alpha - j.beta * j.beta;
    return (a > 0 && b > 0) ? -j.delta * (Math.sqrt(a) - Math.sqrt(b)) : NaN;
  }

  // ---- quatrieme cumulant, par famille ----------------------------------
  // Ne sert qu'a dimensionner la grille de COS et de PROJ.
  function c4Sauts(p, T) {
    let c4 = T * composantes(p).reduce(
      (a, j) => a + j.lam * (j.muJ ** 4 + 6 * j.sigJ ** 2 * j.muJ ** 2 + 3 * j.sigJ ** 4), 0);
    if (p.cgmy && p.cgmy.C > 0) {
      const j = p.cgmy;
      c4 += T * j.C * gammaR(4 - j.Y) * (Math.pow(j.M, j.Y - 4) + Math.pow(j.G, j.Y - 4));
    }
    if (p.kou && p.kou.lam > 0) {
      const j = p.kou;
      c4 += T * j.lam * 24 * (j.p / Math.pow(j.eta1, 4) + (1 - j.p) / Math.pow(j.eta2, 4));
    }
    if (p.vg && p.vg.nu > 0) {
      const j = p.vg, s2 = j.sigma * j.sigma, t2 = j.theta * j.theta;
      c4 += T * (3 * s2 * s2 * j.nu + 12 * s2 * t2 * j.nu * j.nu
                 + 6 * t2 * t2 * Math.pow(j.nu, 3));
    }
    if (p.nig && p.nig.delta > 0) {
      const j = p.nig, g = j.alpha * j.alpha - j.beta * j.beta;
      if (g > 0) c4 += T * 3 * j.delta * j.alpha * j.alpha
                       * (j.alpha * j.alpha + 4 * j.beta * j.beta) / Math.pow(g, 3.5);
    }
    return c4;
  }

'''
anc = "  // log of the Bates futures CF of the log-return (F0 = 1), at real u"
assert s.count(anc) == 1
s = s.replace(anc, NOUVEAU + anc)

# ------------------------------------------------- 2. les branches du logCF
anc = """    if (p.cgmy && p.cgmy.C > 0) {
      const ps = psiCGMY(u, p.cgmy), om = compensateurCGMY(p.cgmy);
      A0 += T * ps[0]; A1 += T * (ps[1] - u * om);
    }"""
assert s.count(anc) == 1
s = s.replace(anc, anc + """
    if (p.kou && p.kou.lam > 0) {
      const ps = psiKou(u, p.kou), om = compensateurKou(p.kou);
      A0 += T * ps[0]; A1 += T * (ps[1] - u * om);
    }
    if (p.vg && p.vg.nu > 0 && p.vg.sigma > 0) {
      const ps = psiVG(u, p.vg), om = compensateurVG(p.vg);
      A0 += T * ps[0]; A1 += T * (ps[1] - u * om);
    }
    if (p.nig && p.nig.delta > 0) {
      const ps = psiNIG(u, p.nig), om = compensateurNIG(p.nig);
      A0 += T * ps[0]; A1 += T * (ps[1] - u * om);
    }""")

# ----------------------------------------- 3. c4 : une seule source de verite
vieux = """    let c4 = T * composantes(p).reduce((a, j) => a + j.lam * (j.muJ ** 4 + 6 * j.sigJ ** 2 * j.muJ ** 2 + 3 * j.sigJ ** 4), 0);
    if (p.cgmy && p.cgmy.C > 0) {
      const j = p.cgmy;
      c4 += T * j.C * gammaR(4 - j.Y) * (Math.pow(j.M, j.Y - 4) + Math.pow(j.G, j.Y - 4));
    }"""
n = s.count(vieux)
assert n >= 1, "bloc c4 introuvable"
s = s.replace(vieux, "    const c4 = c4Sauts(p, T);")

# --------------------------- 3bis. plus de composante lognormale fantome
# `composantes` fabriquait {lam: undefined} des qu'il n'y avait ni sauts ni
# cgmy. Le CF la sautait (`if (!j.lam) continue`), mais le quatrieme cumulant
# la multipliait : c4 devenait NaN et PROJ rendait NaN pour Heston, Kou, VG et
# NIG. Une famille absente doit rendre une liste VIDE, pas une composante nulle.
vieux = ("  const composantes = p => p.sauts || "
         "(p.cgmy ? [] : [{lam: p.lam, muJ: p.muJ, sigJ: p.sigJ}]);")
assert s.count(vieux) == 1, "composantes"
s = s.replace(vieux,
              "  const composantes = p => p.sauts\n"
              "    || (isFinite(p.lam) && p.lam > 0\n"
              "        ? [{lam: p.lam, muJ: p.muJ, sigJ: p.sigJ}] : []);")

# ------------------------------------------------------------ 4. l'interface
anc = "  const api = { logCF, cosMaturity, projMaturity, b76, ivB76, gammaR, psiCGMY, compensateurCGMY };"
assert s.count(anc) == 1
s = s.replace(anc,
              "  const api = { logCF, cosMaturity, projMaturity, b76, ivB76, gammaR,\n"
              "                psiCGMY, compensateurCGMY, psiKou, compensateurKou,\n"
              "                psiVG, compensateurVG, psiNIG, compensateurNIG, c4Sauts };")
s = s.replace("else root.BatesPricer = api;", "else root.Pricer = root.BatesPricer = api;")

io.open(DST / "pricer.js", "w", encoding="utf-8").write(s)
print("pricer unifie ecrit : {} lignes, {} bloc(s) c4 unifie(s)".format(
    s.count("\n") + 1, n))
