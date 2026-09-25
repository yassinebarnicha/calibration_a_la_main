# -*- coding: utf-8 -*-
u"""Ajoute les blocs de VARIANCE au pricer unifie.

Le cadre : log Phi(u) = A_variance + A_sauts_spot, les blocs s'additionnant
dans l'exposant parce que les deux sources sont independantes. Un bloc absent
ne laisse aucun parametre derriere lui -- Black-Scholes n'est pas un Heston a
sigma nul, c'est un exposant qui ne contient pas kappa.

  p.var = {type: "const", sig}                     volatilite constante
        | {type: "heston", kappa, theta, sigma, rho, v0}
        | {type: "heston2", f: [ {..}, {..} ]}     deux facteurs CIR
        | {type: "hestonJ", ..., lamV, muV}        + sauts de variance exponentiels

POURQUOI LES SAUTS DE VARIANCE PASSENT PAR UNE QUADRATURE. Ils n'entrent pas en
addition : ils modifient l'ODE de Riccati, qui gagne le terme

    lamV * integrale_0^T [ E(e^{D(t) Z}) - 1 ] dt ,   Z ~ Exp(muV)

et pour Z exponentielle E(e^{DZ}) = 1/(1 - muV D) est une forme FERMEE. Seule
l'integrale en t reste, sur un integrand lisse : Simpson a 64 points la rend a
la precision machine, pour un cout negligeable. On prefere cette quadrature a
la primitive de Duffie-Pan-Singleton, qu'il faudrait transcrire sans erreur et
qui ne gagnerait rien de mesurable.
"""
import io
from pathlib import Path

P = Path(r"C:/Users/busta/stage/Bates_calibrator/calibration_a_la_main/sources/pricer.js")
s = io.open(P, encoding="utf-8").read()

# ---- 1. on isole le bloc Heston dans une fonction reutilisable -------------
vieux = s[s.index("    const k = p.kappa, th = p.theta"):s.index("    let A0 = 0, A1 = 0;")]
assert vieux.strip().startswith("const k = p.kappa"), "bloc Heston"

NOUVEAU = u'''  // ---- un facteur CIR : le C + D v0 de Heston, isole pour etre reutilise --
  // D(t) est rendu aussi : les sauts de variance en ont besoin le long de [0,T].
  function blocHeston(u, T, q) {
    const k = q.kappa, th = q.theta, s = q.sigma, rho = q.rho;
    const bm = [k, -rho * s * u];
    const qq = cmul([-k, rho * s * u], [-k, rho * s * u]);
    const d = csqrt([qq[0] + s * s * u * u, qq[1] + s * s * u]);
    const num = [bm[0] - d[0], bm[1] - d[1]];
    const den = [bm[0] + d[0], bm[1] + d[1]];
    const g = cdiv(num, den);
    const De = t => {                       // coefficient de v a l'instant t
      const e = cexp([-d[0] * t, -d[1] * t]);
      const ge = cmul(g, e);
      return cmul([num[0] / (s * s), num[1] / (s * s)],
                  cdiv([1 - e[0], -e[1]], [1 - ge[0], -ge[1]]));
    };
    const ed = cexp([-d[0] * T, -d[1] * T]);
    const ged = cmul(g, ed);
    const lg = clog(cdiv([1 - ged[0], -ged[1]], [1 - g[0], -g[1]]));
    const C = [k * th / (s * s) * (num[0] * T - 2 * lg[0]),
               k * th / (s * s) * (num[1] * T - 2 * lg[1])];
    const D = De(T);
    return {C, D, De};
  }

  // ---- sauts de VARIANCE, tailles exponentielles de moyenne muV ----------
  // lamV * int_0^T [ 1/(1 - muV D(t)) - 1 ] dt, Simpson a 64 intervalles.
  function blocSautsVar(T, lamV, muV, De) {
    const n = 64, h = T / n;
    let s0 = 0, s1 = 0;
    for (let i = 0; i <= n; i++) {
      const w = (i === 0 || i === n) ? 1 : (i % 2 ? 4 : 2);
      const D = De(i * h);
      const z = cdiv([1, 0], [1 - muV * D[0], -muV * D[1]]);
      s0 += w * (z[0] - 1); s1 += w * z[1];
    }
    const f = lamV * h / 3;
    return [f * s0, f * s1];
  }

'''

s = s.replace("  // log of the Bates futures CF of the log-return (F0 = 1), at real u",
              NOUVEAU + "  // log of the futures CF of the log-return (F0 = 1), at real u")

# ---- 2. logCF assemble les blocs choisis ----------------------------------
i = s.index("    const k = p.kappa, th = p.theta")
j = s.index("    let A0 = 0, A1 = 0;")
s = s[:i] + '''    // ---- bloc VARIANCE : present seulement s'il est choisi ----
    let V0 = 0, V1 = 0;
    const V = p.var || {type: "heston", kappa: p.kappa, theta: p.theta,
                        sigma: p.sigma, rho: p.rho, v0: p.v0};
    if (V.type === "const") {
      // Black-Scholes : -sig^2 T (u^2 + iu) / 2. Aucun kappa nulle part.
      const s2 = V.sig * V.sig * T / 2;
      V0 = -s2 * u * u; V1 = -s2 * u;
    } else if (V.type === "heston2") {
      // deux facteurs CIR independants : les exposants s'additionnent
      for (const q of V.f) {
        const b = blocHeston(u, T, q);
        V0 += b.C[0] + b.D[0] * q.v0; V1 += b.C[1] + b.D[1] * q.v0;
      }
    } else {
      const b = blocHeston(u, T, V);
      V0 = b.C[0] + b.D[0] * V.v0; V1 = b.C[1] + b.D[1] * V.v0;
      if (V.type === "hestonJ" && V.lamV > 0) {
        const j2 = blocSautsVar(T, V.lamV, V.muV, b.De);
        V0 += j2[0]; V1 += j2[1];
      }
    }
''' + s[j:]

s = s.replace("    return [C[0] + D[0] * v0 + A0, C[1] + D[1] * v0 + A1];",
              "    return [V0 + A0, V1 + A1];")

# c2 de la grille : il vient deja de differences finies sur logCF, rien a faire
s = s.replace("  const api = { logCF,", "  const api = { logCF, blocHeston, blocSautsVar,")

io.open(P, "w", encoding="utf-8").write(s)
print("blocs de variance ajoutes : {} lignes".format(s.count("\n") + 1))
