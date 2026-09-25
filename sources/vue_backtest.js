/* Onglet Backtest : reglages et tableau de bord, mise en page de la methodo.

   Le tableau de bord reprend la structure de la section 4.2 : un tableau de
   synthese par facteur de risque, horizons en colonnes, puis quatre rangees de
   vignettes -- feux tricolores, vue historique, quantiles observes, et ecart
   au median du modele.

   L'onglet est branche SANS toucher au classeur : on observe l'attribut
   aria-selected de la barre d'onglets et on echange la zone centrale quand
   « Backtest » devient actif. Le classeur reste ignorant de son existence. */
(function () {
  "use strict";
  const $ = id => document.getElementById(id);
  const M = window.MoteurBacktest, Q = window.Quantiles, B = window.Backtest;
  const COUL = {vert: "#1a7f4b", ambre: "#d99100", rouge: "#c02040", "—": "#8a8a90"};
  const FOND = {vert: "#e7f5ed", ambre: "#fdf1dc", rouge: "#fae7ec", "—": "#f2f2f4"};

  let vue = null, etat = null, enCours = false;

  // ---------------------------------------------------------------- reglages
  function construireReglages(hote) {
    hote.innerHTML = "";
    const sec = document.createElement("div");
    sec.className = "sec";
    const champ = (lab, html) => {
      const l = document.createElement("label");
      l.className = "field";
      l.innerHTML = "<span>" + lab + "</span>" + html;
      sec.append(l);
      return l;
    };
    const dates = JOURS.map(j => j.date);
    const opts = (sel) => dates.map(d =>
      '<option value="' + d + '"' + (d === sel ? " selected" : "") + ">" + fdate(d) + "</option>").join("");
    const i60 = dates[Math.max(0, dates.length - 60)], dern = dates[dates.length - 1];

    champ("Paramètres — de", '<select id="bt-p0">' + opts(i60) + "</select>");
    champ("Paramètres — à", '<select id="bt-p1">' + opts(dern) + "</select>");
    champ("Quels paramètres",
      '<select id="bt-src">'
      + '<option value="ecran">Ceux de l’écran, figés</option>'
      + '<option value="jour">Le réglage enregistré du jour</option>'
      + "</select>");
    champ("Si plusieurs par jour",
      '<select id="bt-choix">'
      + '<option value="rmse">le meilleur RMSE IV</option>'
      + '<option value="dernier">le plus récent</option>'
      + '<option value="final">celui marqué ★ final</option>'
      + "</select>");
    champ("… filtrer sur le nombre d’échéances",
      '<input class="num" id="bt-nech" type="number" min="0" step="1" placeholder="toutes">');
    champ("… filtrer sur la perte",
      '<select id="bt-perte"><option value="">toutes</option>'
      + Object.keys(MESURES).map(k => '<option value="' + k + '">' + MESURES[k].court + "</option>").join("")
      + "</select>");

    const sep = document.createElement("div");
    sep.className = "hint";
    sep.style.borderTop = "1px solid var(--line)";
    sep.style.paddingTop = "8px";
    sep.innerHTML = "La plage des <b>paramètres</b> dit d’où viennent les jeux calibrés. "
      + "La plage du <b>backtesting</b> dit sur quelle fenêtre on observe. Elles sont "
      + "indépendantes : on peut calibrer sur août et observer septembre.";
    sec.append(sep);

    champ("Backtesting — de", '<select id="bt-b0">' + opts(i60) + "</select>");
    champ("Backtesting — à", '<select id="bt-b1">' + opts(dern) + "</select>");
    champ("Horizons (jours, séparés par des virgules)",
      '<input class="num" id="bt-h" type="text" value="1, 5, 10, 21">');
    champ("Pas entre les ancres (jours)",
      '<input class="num" id="bt-pas" type="number" min="1" step="1" value="1">');
    champ("Tolérance sur la date cible (jours)",
      '<input class="num" id="bt-tol" type="number" min="0" max="5" step="1" value="1">');
    champ("Intervalles des feux",
      '<input class="num" id="bt-int" type="number" min="2" max="50" step="1" value="10">');

    const lf = document.createElement("div");
    lf.className = "field";
    lf.innerHTML = "<span>Facteurs de risque</span>"
      + '<div id="bt-fact" style="display:grid;gap:3px;max-height:190px;overflow:auto;'
      + 'border:1px solid var(--line);padding:6px"></div>';
    sec.append(lf);

    const row = document.createElement("div");
    row.className = "row";
    row.innerHTML = '<button class="btn primary" id="bt-go" type="button">Lancer le backtest</button>'
      + '<button class="btn" id="bt-csv" type="button">Exporter les rangs</button>';
    sec.append(row);
    const st = document.createElement("div");
    st.className = "status"; st.id = "bt-etat";
    st.textContent = "Choisis les plages, puis lance.";
    sec.append(st);
    hote.append(sec);

    majFacteurs();
    $("bt-b0").addEventListener("change", majFacteurs);
    $("bt-go").addEventListener("click", lancer);
    $("bt-csv").addEventListener("click", exporter);
  }

  /* Les facteurs proposes sont les echeances cotees au debut de la fenetre
     d'observation, plus le comptant. Une echeance qui expire pendant la
     fenetre reste proposee : le moteur s'arrete simplement quand elle
     disparait des instantanes. */
  function majFacteurs() {
    const z = $("bt-fact");
    if (!z) return;
    const d0 = $("bt-b0").value;
    const src = CACHE[d0] || D;
    const ech = [...new Set(src.mats.map(m => m.date))].sort();
    z.innerHTML = '<label class="chk"><input type="checkbox" id="bt-f-spot" checked> '
      + "<b>comptant</b> (déduit de r)</label>"
      + ech.map((e, i) => '<label class="chk"><input type="checkbox" class="bt-f-ech" '
          + 'value="' + e + '"' + (i >= ech.length - 4 ? " checked" : "") + "> future " + fdate(e)
          + "</label>").join("");
  }

  // ------------------------------------------------------------- lancement
  function jeuxDuJour(jour) {
    const loc = (typeof lireReg === "function" ? (lireReg()[jour] || []) : []);
    const dsk = (typeof DISQUE !== "undefined" ? (DISQUE[jour] || []) : []);
    return [...loc, ...dsk];
  }

  function choisirJeu(jour) {
    const nech = parseInt($("bt-nech").value, 10);
    const perte = $("bt-perte").value;
    let l = jeuxDuJour(jour).filter(r => r && r.parametres);
    if (isFinite(nech)) l = l.filter(r => (r.donnees && r.donnees.echeances || []).length === nech);
    if (perte) l = l.filter(r => r.perte && r.perte.mesure === perte);
    if (!l.length) return null;
    const mode = $("bt-choix").value;
    if (mode === "final") { const f = l.find(r => r.final); return f ? f.parametres : null; }
    if (mode === "dernier") return l[l.length - 1].parametres;
    const bps = r => (r.resultats && r.resultats.rmse_iv_bps != null) ? r.resultats.rmse_iv_bps : 1e9;
    return l.reduce((a, b2) => bps(b2) < bps(a) ? b2 : a).parametres;
  }

  async function lancer() {
    if (enCours) return;
    enCours = true;
    const dire = t => { const e = $("bt-etat"); if (e) e.textContent = t; };
    try {
      const p0 = $("bt-p0").value, p1 = $("bt-p1").value;
      const b0 = $("bt-b0").value, b1 = $("bt-b1").value;
      const pas = Math.max(1, +$("bt-pas").value || 1);
      const tol = Math.max(0, +$("bt-tol").value || 0);
      const inter = Math.max(2, +$("bt-int").value || 10);
      const horizons = $("bt-h").value.split(",").map(x => parseInt(x, 10))
        .filter(x => isFinite(x) && x > 0);
      if (!horizons.length) { dire("Aucun horizon valide."); return; }

      const dates = JOURS.map(j => j.date);
      // une ancre doit etre dans LES DEUX plages : elle fournit les parametres
      // et ouvre une fenetre d'observation
      let ancres = dates.filter(d => d >= p0 && d <= p1 && d >= b0 && d <= b1);
      ancres = ancres.filter((d, i) => i % pas === 0);
      if (!ancres.length) { dire("Les deux plages ne se croisent sur aucun jour."); return; }

      const facteurs = [];
      if ($("bt-f-spot").checked) facteurs.push({type: "spot"});
      document.querySelectorAll(".bt-f-ech:checked").forEach(c =>
        facteurs.push({type: "future", ech: c.value}));
      if (!facteurs.length) { dire("Choisis au moins un facteur de risque."); return; }

      const fige = $("bt-src").value === "ecran";
      const courant = {var: BV().faire(S.p), ...BS().faire(S.p)};
      const params = jour => {
        if (fige) return courant;
        const p = choisirJeu(jour);
        if (!p) return null;
        // un jeu enregistre porte les noms de SA composition : on ne le
        // convertit pas, on l'assemble avec la composition courante
        return {var: BV().faire(p), ...BS().faire(p)};
      };

      dire("Chargement… 0 / " + ancres.length);
      const res = await M.lancer({
        ancres, params, horizons, facteurs, tolerance: tol,
        disponibles: dates.filter(d => d >= b0 && d <= b1),
        charger: async j => CACHE[j] || await chargerFichier(j),
      }, (k, n) => dire("Calcul… " + k + " / " + n + " ancres"));

      etat = {res, horizons, pas, inter, ancres, fige};
      dessiner();
      const n = res.detail.length;
      dire(n + " observation(s) sur " + ancres.length + " ancres, "
           + facteurs.length + " facteur(s), " + horizons.length + " horizon(s).");
    } catch (e) {
      dire("Échec : " + e.message);
    } finally { enCours = false; }
  }

  function exporter() {
    if (!etat) return;
    const L = ["facteur;ancre;cible;horizon;horizon_reel;base;observe;log_rendement;rang"];
    etat.res.detail.forEach(d => L.push([d.facteur, d.ancre, d.cible, d.horizon,
      d.horizonReel, d.base, d.obs, d.logRendement.toFixed(8), d.rang.toFixed(6)].join(";")));
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([L.join("\n")], {type: "text/csv"}));
    a.download = "backtest_rangs_" + new Date().toISOString().slice(0, 10) + ".csv";
    document.body.appendChild(a); a.click(); a.remove();
  }

  // ----------------------------------------------------------- tableau de bord
  function dessiner() {
    if (!vue || !etat) return;
    vue.innerHTML = "";
    const {res, horizons, pas, inter} = etat;
    res.cles.forEach(cle => {
      const bloc = document.createElement("div");
      bloc.className = "panel";
      bloc.style.padding = "12px 14px";
      const h = document.createElement("h3");
      h.style.margin = "0 0 8px";
      h.textContent = "Facteur de risque : " + cle;
      bloc.append(h);

      const anas = {};
      horizons.forEach(hz => {
        const r = res.par[cle][hz] || [];
        anas[hz] = r.length >= 3
          ? B.analyser(r, {intervalles: inter, recouvrement: hz > pas, lags: M.lagsPour(hz, pas)})
          : null;
      });
      bloc.append(tableau(horizons, anas));
      bloc.append(vignettes(cle, horizons, anas));
      vue.append(bloc);
    });
  }

  function cellule(txt, zone) {
    const td = document.createElement("td");
    td.textContent = txt;
    td.style.padding = "3px 8px";
    td.style.textAlign = "right";
    td.style.fontFamily = "var(--mono)";
    td.style.fontSize = "12px";
    if (zone) { td.style.background = FOND[zone]; td.style.color = COUL[zone]; }
    return td;
  }

  function tableau(horizons, anas) {
    const t = document.createElement("table");
    t.style.borderCollapse = "collapse";
    t.style.width = "100%";
    t.style.marginBottom = "10px";
    const entete = t.insertRow();
    const th0 = document.createElement("th");
    th0.style.textAlign = "left"; th0.style.padding = "3px 8px"; th0.style.fontSize = "12px";
    entete.append(th0);
    horizons.forEach(h => {
      const th = document.createElement("th");
      th.textContent = h + (h > 1 ? " jours" : " jour");
      th.style.padding = "3px 8px"; th.style.fontSize = "12px";
      entete.append(th);
    });
    const ligne = (lab, f) => {
      const tr = t.insertRow();
      const td = document.createElement("td");
      td.textContent = lab;
      td.style.padding = "3px 8px"; td.style.fontSize = "12px"; td.style.whiteSpace = "nowrap";
      tr.append(td);
      horizons.forEach(h => tr.append(f(anas[h])));
    };
    const pc = (v, d) => isFinite(v) ? (v * 100).toFixed(d == null ? 1 : d) + " %" : "—";
    ligne("Nombre de points", a => cellule(a ? a.n : "—"));
    ligne("Recouvrement / dépendance", a => a && a.ljungBox
      ? cellule(pc(a.ljungBox.p, 2), a.ljungBox.zone) : cellule("aucun"));
    ligne("Test d’adéquation : KS", a => a ? cellule(pc(a.ks.p, 2), a.ks.zone) : cellule("—"));
    ligne("Test d’adéquation : AD", a => a ? cellule(pc(a.ad.p, 2), a.ad.zone) : cellule("—"));
    ligne("Moyenne : observée / attendue", a => a ? cellule(pc(a.moyenneRatio)) : cellule("—"));
    ligne("Test d’égalité des médianes", a => a ? cellule(pc(a.signe.p, 2), a.signe.zone) : cellule("—"));
    ligne("Variance : observée / attendue", a => a ? cellule(pc(a.variance.ratio)) : cellule("—"));
    ligne("Test d’égalité des variances", a => a ? cellule(pc(a.variance.p, 2), a.variance.zone) : cellule("—"));
    ligne("Feux tricolores", a => a ? cellule(a.feux.zone, a.feux.zone) : cellule("—"));
    return t;
  }

  function vignettes(cle, horizons, anas) {
    const g = document.createElement("div");
    g.style.display = "grid";
    g.style.gridTemplateColumns = "150px repeat(" + horizons.length + ", minmax(0,1fr))";
    g.style.gap = "6px";
    g.style.alignItems = "center";
    const RANGEES = [
      ["Validation<br>feux tricolores", (a, cv) => traceFeux(cv, a)],
      ["Vue historique<br>réalisé vs bande", (a, cv, h) => traceHisto(cv, cle, h)],
      ["Quantiles observés", (a, cv) => traceRangs(cv, a)],
      ["Réalisé moins<br>médiane du modèle", (a, cv, h) => traceEcart(cv, cle, h)],
    ];
    RANGEES.forEach(([lab, f]) => {
      const t = document.createElement("div");
      t.innerHTML = lab;
      t.style.fontSize = "11px"; t.style.color = "var(--ink2)"; t.style.textAlign = "right";
      g.append(t);
      horizons.forEach(h => {
        const box = document.createElement("div");
        box.style.height = "92px";
        box.style.border = "1px solid var(--line)";
        box.style.background = "var(--panel)";
        const cv = document.createElement("canvas");
        cv.style.width = "100%"; cv.style.height = "100%"; cv.style.display = "block";
        box.append(cv); g.append(box);
        requestAnimationFrame(() => f(anas[h], cv, h));
      });
    });
    return g;
  }

  function ctx(cv) {
    const dpr = window.devicePixelRatio || 1;
    const W = cv.clientWidth, H = cv.clientHeight;
    if (!W || !H) return null;
    cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
    const g = cv.getContext("2d");
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, W, H);
    return {g, W, H};
  }

  /* Feux : bandes RAG en fond, ratio realise/attendu par intervalle par-dessus.
     Les bandes reprennent l'echelle de la methodo, 0 a 2.5. */
  function traceFeux(cv, a) {
    const c = ctx(cv); if (!c || !a) return;
    const {g, W, H} = c, HAUT = 2.5;
    const y = v => H - (v / HAUT) * H;
    // vert au centre, ambre autour, rouge aux extremes
    [[0, 0.5, "#c02040"], [0.5, 0.75, "#d99100"], [0.75, 1.35, "#1a7f4b"],
     [1.35, 2.0, "#d99100"], [2.0, HAUT, "#c02040"]].forEach(([lo, hi, col]) => {
      g.fillStyle = col; g.globalAlpha = .82;
      g.fillRect(0, y(hi), W, y(lo) - y(hi));
    });
    g.globalAlpha = 1;
    g.strokeStyle = "#111"; g.lineWidth = 1.4;
    g.beginPath();
    a.feux.cases.forEach((x, i) => {
      const px = (i + 0.5) / a.feux.cases.length * W;
      const py = Math.max(2, Math.min(H - 2, y(isFinite(x.ratio) ? x.ratio : 0)));
      i ? g.lineTo(px, py) : g.moveTo(px, py);
    });
    g.stroke();
    g.fillStyle = "#111";
    a.feux.cases.forEach((x, i) => {
      const px = (i + 0.5) / a.feux.cases.length * W;
      const py = Math.max(2, Math.min(H - 2, y(isFinite(x.ratio) ? x.ratio : 0)));
      g.beginPath(); g.arc(px, py, 1.8, 0, 2 * Math.PI); g.fill();
    });
  }

  /* Vue historique : le log-rendement realise, avec la bande des quantiles. */
  function traceHisto(cv, cle, h) {
    const c = ctx(cv); if (!c || !etat) return;
    const {g, W, H} = c;
    const d = etat.res.detail.filter(x => x.facteur === cle && x.horizon === h);
    if (d.length < 2) return;
    const vals = d.map(x => x.logRendement);
    let lo = Math.min(...vals), hi = Math.max(...vals);
    const m = (hi - lo) * .2 || .01; lo -= m; hi += m;
    const y = v => H - (v - lo) / (hi - lo) * H;
    g.strokeStyle = "#999"; g.setLineDash([2, 2]); g.lineWidth = 1;
    g.beginPath(); g.moveTo(0, y(0)); g.lineTo(W, y(0)); g.stroke();
    g.setLineDash([]);
    g.strokeStyle = "#3A0850"; g.lineWidth = 1.6;
    g.beginPath();
    d.forEach((x, i) => {
      const px = d.length > 1 ? i / (d.length - 1) * W : W / 2;
      i ? g.lineTo(px, y(x.logRendement)) : g.moveTo(px, y(x.logRendement));
    });
    g.stroke();
  }

  /* Quantiles observes : le rang au fil du temps, entre 0 et 1. */
  function traceRangs(cv, a) {
    const c = ctx(cv); if (!c || !a) return;
    const {g, W, H} = c;
    g.strokeStyle = "#ddd"; g.lineWidth = 1;
    [0.25, 0.5, 0.75].forEach(q => {
      g.beginPath(); g.moveTo(0, H * (1 - q)); g.lineTo(W, H * (1 - q)); g.stroke();
    });
    g.strokeStyle = "#c02040"; g.lineWidth = 1.2;
    g.beginPath();
    a.rangs.forEach((r, i) => {
      const px = a.rangs.length > 1 ? i / (a.rangs.length - 1) * W : W / 2;
      const py = H * (1 - r);
      i ? g.lineTo(px, py) : g.moveTo(px, py);
    });
    g.stroke();
  }

  /* Ecart au median du modele, en batonnets. */
  function traceEcart(cv, cle, h) {
    const c = ctx(cv); if (!c || !etat) return;
    const {g, W, H} = c;
    const d = etat.res.detail.filter(x => x.facteur === cle && x.horizon === h);
    if (!d.length) return;
    // la mediane du modele est le quantile 50 % : on la relit du rang
    const ec = d.map(x => x.rang - 0.5);
    const mx = Math.max(...ec.map(Math.abs)) || 0.5;
    const y0 = H / 2;
    g.strokeStyle = "#999"; g.lineWidth = 1;
    g.beginPath(); g.moveTo(0, y0); g.lineTo(W, y0); g.stroke();
    g.strokeStyle = "#3A0850"; g.lineWidth = 1;
    ec.forEach((v, i) => {
      const px = ec.length > 1 ? (i + 0.5) / ec.length * W : W / 2;
      g.beginPath(); g.moveTo(px, y0); g.lineTo(px, y0 - v / mx * (H / 2 - 3)); g.stroke();
    });
  }

  // --------------------------------------------------------------- branchement
  function brancher() {
    const nav = document.querySelector(".onglets");
    const espace = document.querySelector(".espace");
    const principal = espace && espace.querySelector("main");
    if (!nav || !espace || !principal) return;

    vue = document.createElement("div");
    vue.className = "vue-backtest";
    vue.hidden = true;
    vue.style.display = "grid";
    vue.style.gap = "12px";
    vue.style.alignContent = "start";
    vue.innerHTML = '<div class="panel" style="padding:14px 16px">'
      + "<b>Backtest du facteur de risque.</b> Règle les deux plages à droite, "
      + "puis lance. Les rangs sont calculés sur 101 quantiles archivés, obtenus "
      + "par inversion de Fourier — pas de Monte-Carlo, donc pas de rang saturé "
      + "et deux passages donnent le même résultat.</div>";
    espace.insertBefore(vue, espace.querySelector("aside.droite"));
    // la vue occupe la meme case que les smiles
    vue.style.gridColumn = "1"; vue.style.gridRow = "1";

    const panneau = document.querySelector('#bt-hote');
    if (panneau) construireReglages(panneau);

    // on observe l'onglet actif plutot que de modifier le classeur
    const boutons = [...nav.querySelectorAll(".onglet")];
    const iBt = boutons.findIndex(b => /backtest/i.test(b.textContent));
    if (iBt < 0) return;
    const appliquer = () => {
      const actif = boutons[iBt].getAttribute("aria-selected") === "true";
      vue.hidden = !actif;
      const autreVue = document.querySelector(".vue-futures");
      if (actif) {
        principal.hidden = true;
        if (autreVue) autreVue.hidden = true;
        if (etat) dessiner();
      } else if (!autreVue || autreVue.hidden) {
        principal.hidden = false;
      }
    };
    new MutationObserver(appliquer).observe(nav, {subtree: true, attributes: true,
      attributeFilter: ["aria-selected"]});
    appliquer();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", brancher);
  else setTimeout(brancher, 0);
})();
