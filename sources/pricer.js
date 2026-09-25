/* Pricer unifie : Heston + une famille de sauts au choix, COS et PROJ.

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
(function (root) {
  "use strict";

  // ---- complex helpers on [re, im] pairs --------------------------------
  function cmul(a, b) { return [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]; }
  function cdiv(a, b) {
    const d = b[0] * b[0] + b[1] * b[1];
    return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d];
  }
  function cexp(a) { const e = Math.exp(a[0]); return [e * Math.cos(a[1]), e * Math.sin(a[1])]; }
  function clog(a) { return [Math.log(Math.hypot(a[0], a[1])), Math.atan2(a[1], a[0])]; }
  function csqrt(a) {                       // principal branch, like numpy
    const r = Math.hypot(a[0], a[1]);
    const re = Math.sqrt((r + a[0]) / 2);
    let im = Math.sqrt(Math.max((r - a[0]) / 2, 0));
    if (a[1] < 0) im = -im;
    return [re, im];
  }

  const composantes = p => p.sauts
    || (isFinite(p.lam) && p.lam > 0
        ? [{lam: p.lam, muJ: p.muJ, sigJ: p.sigJ}] : []);

  // z^y pour z complexe, y reel : exp(y log z), branche principale
  function cpow(z, y) { const l = clog(z); return cexp([y * l[0], y * l[1]]); }

  // Gamma d'Euler, Lanczos g=7 n=9, puis reflexion pour les arguments negatifs.
  // CGMY a besoin de Gamma(-Y) avec 0 < Y < 2 : toujours negatif non entier.
  const LANCZOS = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                   771.32342877765313, -176.61502916214059, 12.507343278686905,
                   -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  function gammaR(x) {
    if (x < 0.5) return Math.PI / (Math.sin(Math.PI * x) * gammaR(1 - x));
    x -= 1;
    let a = LANCZOS[0];
    const t = x + 7.5;
    for (let i = 1; i < 9; i++) a += LANCZOS[i] / (x + i);
    return Math.sqrt(2 * Math.PI) * Math.pow(t, x + 0.5) * Math.exp(-t) * a;
  }

  // Exposant de saut CGMY : psi(u) = C Gamma(-Y) [(M-iu)^Y - M^Y + (G+iu)^Y - G^Y].
  // Le compensateur impose A(-i) = 0, donc E[F_T] = F_0 : on retranche iu psi(-i),
  // et psi(-i) est REEL car (M - i(-i)) = M-1 et (G + i(-i)) = G+1.
  // Gamma(-Y) a des poles en Y = 0 et Y = 1 : on s'en ecarte d'un epsilon
  // plutot que de rendre NaN en silence. Y = 0, c'est Variance Gamma, qui a
  // sa propre entree ; Y = 1 n'a pas de sens limite utile.
  const YsansPole = y => (Math.abs(y) < 1e-4 ? (y < 0 ? -1e-4 : 1e-4)
    : Math.abs(y - 1) < 1e-4 ? (y < 1 ? 1 - 1e-4 : 1 + 1e-4) : y);
  function psiCGMY(u, j) {
    const Y = YsansPole(j.Y), gY = gammaR(-Y);
    const a = cpow([j.M, -u], Y);            // (M - iu)^Y
    const b = cpow([j.G, u], Y);             // (G + iu)^Y
    const cst = Math.pow(j.M, Y) + Math.pow(j.G, Y);
    return [j.C * gY * (a[0] + b[0] - cst), j.C * gY * (a[1] + b[1])];
  }
  function compensateurCGMY(j) {
    const Y = YsansPole(j.Y), gY = gammaR(-Y);
    return j.C * gY * (Math.pow(j.M - 1, Y) - Math.pow(j.M, Y)
                     + Math.pow(j.G + 1, Y) - Math.pow(j.G, Y));
  }


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
      const Yc = YsansPole(j.Y);
      c4 += T * j.C * gammaR(4 - Yc) * (Math.pow(j.M, Yc - 4) + Math.pow(j.G, Yc - 4));
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

  // ---- un facteur CIR : le C + D v0 de Heston, isole pour etre reutilise --
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

  // log of the futures CF of the log-return (F0 = 1), at real u
  function logCF(u, T, p) {
    // ---- bloc VARIANCE : present seulement s'il est choisi ----
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
    let A0 = 0, A1 = 0;
    for (const j of composantes(p)) {
      if (!j.lam) continue;
      const kJ = Math.exp(j.muJ + 0.5 * j.sigJ * j.sigJ) - 1;
      const ej = cexp([-0.5 * j.sigJ * j.sigJ * u * u, u * j.muJ]);
      A0 += j.lam * T * (ej[0] - 1); A1 += j.lam * T * (ej[1] - u * kJ);
    }
    if (p.cgmy && p.cgmy.C > 0) {
      const ps = psiCGMY(u, p.cgmy), om = compensateurCGMY(p.cgmy);
      A0 += T * ps[0]; A1 += T * (ps[1] - u * om);
    }
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
    }
    return [V0 + A0, V1 + A1];
  }

  // prices (USD) for all rows of one maturity; rows: {F, K, r, c}
  function cosMaturity(rows, T, p, N, L) {
    N = N || 512; L = L || 24;
    const h = 1e-5;
    const l0 = logCF(0, T, p), lp = logCF(h, T, p), lm = logCF(-h, T, p);
    // c1 = Re[(lp - lm) / (2 i h)] = Im(lp - lm) / (2h) ; c2 = -Re(lp - 2 l0 + lm)/h^2
    const c1 = (lp[1] - lm[1]) / (2 * h);
    const c2 = -(lp[0] - 2 * l0[0] + lm[0]) / (h * h);
    const out = new Float64Array(rows.length);
    if (!isFinite(c1) || !isFinite(c2) || c2 <= 0) { out.fill(NaN); return out; }
    const width = 2 * L * Math.sqrt(c2);
    const shift = c1 - 0.5 * width;
    const w = new Float64Array(N), re = new Float64Array(N);
    for (let k = 0; k < N; k++) {
      const u = k * Math.PI / width;
      w[k] = u;
      const lc = logCF(u, T, p);
      // phi(u) * exp(-i u shift): only the real part is needed (Vk is real)
      const ph = cexp([lc[0], lc[1] - u * shift]);
      re[k] = isFinite(ph[0]) ? ph[0] : NaN;
    }
    for (let j = 0; j < rows.length; j++) {
      const R = rows[j];
      const a = Math.log(R.F / R.K) + shift;
      const ea = Math.exp(a);
      let acc = 0;
      for (let k = 0; k < N; k++) {
        const u = w[k];
        let chi, psi;
        if (k === 0) { chi = 1 - ea; psi = -a; }
        else {
          const cu = Math.cos(u * a), su = Math.sin(u * a);   // cos(-ua)=cu, sin(-ua)=-su
          chi = (cu - u * su - ea) / (1 + u * u);
          psi = -su / u;
        }
        let v = (2 / width) * R.K * (psi - chi);
        if (k === 0) v *= 0.5;
        acc += v * re[k];
      }
      const disc = Math.exp(-R.r * T);
      const put = disc * acc;
      out[j] = R.c ? put + disc * (R.F - R.K) : put;
    }
    return out;
  }

  // ---- Black-76 ----------------------------------------------------------
  function ncdf(x) {                      // erfc de Numerical Recipes, ~1.2e-7 relatif (IV a 0.01 bp pres)
    const z = Math.abs(x) / Math.SQRT2;
    const t = 1 / (1 + 0.5 * z);
    const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
      t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
      t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? 1 - 0.5 * r : 0.5 * r;
  }
  function b76(F, K, T, r, sig, call) {
    const df = Math.exp(-r * T);
    if (sig <= 0) return df * Math.max(call ? F - K : K - F, 0);
    const sq = sig * Math.sqrt(T);
    const d1 = (Math.log(F / K) + 0.5 * sq * sq) / sq, d2 = d1 - sq;
    return call ? df * (F * ncdf(d1) - K * ncdf(d2)) : df * (K * ncdf(-d2) - F * ncdf(-d1));
  }
  function ivB76(price, F, K, T, r, call) {
    const df = Math.exp(-r * T);
    const intr = df * Math.max(call ? F - K : K - F, 0);
    const up = df * (call ? F : K);
    if (!(isFinite(price)) || price < intr - 1e-10 || price > up + 1e-10) return NaN;
    if (Math.abs(price - intr) < 1e-12) return 0;
    let lo = 1e-6, hi = 5;
    if (b76(F, K, T, r, hi, call) < price) return NaN;
    let s = 0.6;
    for (let i = 0; i < 100; i++) {
      const v = b76(F, K, T, r, s, call) - price;
      if (Math.abs(v) < 1e-10) return s;
      if (v > 0) hi = s; else lo = s;
      const sq = Math.sqrt(T), d1 = (Math.log(F / K) + 0.5 * s * s * T) / (s * sq);
      const vega = df * F * Math.exp(-0.5 * d1 * d1) / Math.sqrt(2 * Math.PI) * sq;
      let sn = s - v / vega;
      if (!(sn > lo && sn < hi) || !isFinite(sn)) sn = 0.5 * (lo + hi);
      if (Math.abs(sn - s) < 1e-13) return sn;
      s = sn;
    }
    return s;
  }

  // ---- PROJ (Kirkby 2015, spline lineaire, grille decalee) ------------------
  // Port de pricer_proj_vec._proj_prices_one_maturity : une FFT par echeance,
  // partagee par tous les strikes ; F0 = F de la premiere ligne de l'echeance.
  function fft(re, im) {                     // radix-2, noyau e^{-2 pi i kn/N} (numpy)
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const half = len >> 1, ang = -2 * Math.PI / len;
      for (let j = 0; j < half; j++) {
        const wr = Math.cos(ang * j), wi = Math.sin(ang * j);
        for (let i = j; i < n; i += len) {
          const k = i + half;
          const vr = re[k] * wr - im[k] * wi, vi = re[k] * wi + im[k] * wr;
          re[k] = re[i] - vr; im[k] = im[i] - vi; re[i] += vr; im[i] += vi;
        }
      }
    }
  }
  function projMaturity(rows, T, p, N, L) {
    N = N || 2048; L = L || 20;
    const h = 1e-5, n = rows.length, out = new Float64Array(n);
    const F0 = rows[0].F;
    const l0 = logCF(0, T, p), lp = logCF(h, T, p), lm = logCF(-h, T, p);
    const c1 = (lp[1] - lm[1]) / (2 * h), c2 = -(lp[0] - 2 * l0[0] + lm[0]) / (h * h);
    if (!isFinite(c1) || !isFinite(c2) || c2 <= 0) { out.fill(NaN); return out; }
    // c4 fixe la LARGEUR de la grille : une mauvaise valeur ne fait pas diverger
    // le prix, elle degrade silencieusement sa precision. Formule fermee par famille.
    const c4 = c4Sauts(p, T);
    const lws = rows.map(R => Math.log(R.K / F0));
    let maxAbs = 0, maxL = -Infinity;
    lws.forEach(v => { maxAbs = Math.max(maxAbs, Math.abs(v)); maxL = Math.max(maxL, v); });
    const alph = Math.max(L * Math.sqrt(Math.abs(c2) + Math.sqrt(Math.abs(c4))), 1.15 * maxAbs + c1);
    const dx = 2 * alph / (N - 1), a = 1 / dx, lam0 = c1 - (N / 2 - 1) * dx;
    const nbar = (x, ref) => Math.min(Math.floor(a * (x - ref) + 1), N - 1);
    const maxNbar = nbar(maxL, lam0);
    const e75 = Math.exp(-0.75 * dx), e5 = Math.exp(-0.5 * dx), e25 = Math.exp(-0.25 * dx);
    const g1 = 0.5 - (1 / 15) * (7 / 6 + 4 / 3 * e75 + e5 + 4 * e25);
    const g2 = (7 / 3 + 8 / 3 * Math.cosh(0.75 * dx) + 2 * Math.cosh(0.5 * dx) + 8 * Math.cosh(0.25 * dx)) / 15;
    const th0 = (1 / 15) * (7 / 6 + 4 / 3 * e75 + e5 + 4 * e25);
    const qP = 0.5 * (1 + Math.sqrt(3 / 5)), qM = 0.5 * (1 - Math.sqrt(3 / 5));
    const dw = 2 * Math.PI / (N * dx), cons = 24 * a * a;
    const xmin = maxL - (maxNbar - 1) * dx;
    const re = new Float64Array(N), im = new Float64Array(N);
    re[0] = 1 / cons;
    for (let k = 1; k < N; k++) {
      const w = dw * k, lc = logCF(w, T, p);
      const z = (Math.sin(w / (2 * a)) / w) ** 2 / (2 + Math.cos(w / a));
      const mag = z * Math.exp(lc[0]), ph = lc[1] - xmin * w;
      re[k] = mag * Math.cos(ph); im[k] = mag * Math.sin(ph);
      if (!isFinite(re[k]) || !isFinite(im[k])) { out.fill(NaN); return out; }
    }
    fft(re, im);
    const beta = re;
    const nExp = Math.max(maxNbar - 1, 1);
    // garde-fou du Python : grille trop large pour le float64 -> NaN, jamais 0 silencieux
    if (!isFinite(Math.exp(dx * (nExp - 1)))) { out.fill(NaN); return out; }
    const cumB = new Float64Array(N + 1), cumBE = new Float64Array(N + 1);
    for (let j = 0; j < N; j++) {
      cumB[j + 1] = cumB[j] + beta[j];
      cumBE[j + 1] = cumBE[j] + (j < nExp ? beta[j] * Math.exp(dx * j) : 0);
    }
    const exmin = Math.exp(xmin), cl = i => beta[Math.min(Math.max(i, 0), N - 1)];
    for (let j = 0; j < n; j++) {
      const R = rows[j], W = R.K, nb = nbar(lws[j], xmin), m = Math.max(nb - 1, 0);
      const plate = W * cumB[m] - F0 * exmin * g2 * cumBE[m];
      const rho = lws[j] - (xmin + (nb - 1) * dx), z = a * rho;
      const zP = z * qP, zM = z * qM, rP = rho * qP, rM = rho * qM;
      const db0 = z * (1 - 0.5 * z), dbP1 = z - db0;
      const d0 = (z / 18) * (4 * (2 - z) * Math.exp(rho / 2) + 5 * ((1 - zM) * Math.exp(rM) + (1 - zP) * Math.exp(rP)));
      const dP1 = (z / 18) * Math.exp(-dx) * (4 * z * Math.exp(0.5 * rho) + 5 * (zM * Math.exp(rM) + zP * Math.exp(rP)));
      const Gnb = W * (dbP1 - Math.exp(-rho) * Math.exp(dx) * dP1);
      const Gnb1 = W * g1 + W * (db0 - Math.exp(-rho) * (th0 + d0) + th0);
      const tot = plate + cl(nb - 1) * Gnb1 + cl(nb) * Gnb;
      const disc = Math.exp(-R.r * T);
      let px = (cons / N) * disc * tot;
      if (R.c) px += disc * (F0 - W);
      out[j] = Math.max(0, px);
    }
    return out;
  }

  const api = { logCF, blocHeston, blocSautsVar, cosMaturity, projMaturity, b76, ivB76, gammaR,
                psiCGMY, compensateurCGMY, psiKou, compensateurKou,
                psiVG, compensateurVG, psiNIG, compensateurNIG, c4Sauts };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.Pricer = root.BatesPricer = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
