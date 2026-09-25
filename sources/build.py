# -*- coding: utf-8 -*-
u"""Construit ../index.html : la page, le pricer et un jour integres.

Les autres jours restent dans ../jours/<date>.js et sont charges a la demande
par une balise <script> classique -- pas par fetch, qui serait bloque en
file://. C'est ce qui permet au dossier de marcher par double-clic, sans
serveur, sur un poste isole.

Usage : python sources/build.py
"""
import io
import json
import os
from pathlib import Path

ICI = Path(__file__).resolve().parent
RACINE = ICI.parent
PRE, SUF = "window.__recevoirJour(", ");"
PRE_FUT, SUF_FUT = "window.__recevoirFutures(", ");"


def jour_json(date):
    brut = (RACINE / "jours" / "{}.js".format(date)).read_text(encoding="utf-8")
    assert brut.startswith(PRE) and brut.endswith(SUF), date
    return brut[len(PRE):-len(SUF)]


def futures_json():
    p = RACINE / "futures.js"
    if not p.exists():
        return "null"
    brut = p.read_text(encoding="utf-8")
    assert brut.startswith(PRE_FUT) and brut.endswith(SUF_FUT), "futures.js"
    return brut[len(PRE_FUT):-len(SUF_FUT)]


def index_jours():
    u"""date, nombre d'options et d'echeances : de quoi peupler la liste."""
    out = []
    for p in sorted((RACINE / "jours").glob("*.js")):
        d = json.loads(jour_json(p.stem))
        out.append({"date": d["jour"], "n": len(d["cols"]["m"]),
                    "nprod": sum(1 for x in d["cols"]["prod"] if x),
                    "nmat": len(d["mats"])})
    return out


def main():
    index = index_jours()
    assert index, "aucun jour dans jours/"
    jour0 = index[-1]["date"]          # le plus recent par defaut
    html = (ICI / "template_full.html").read_text(encoding="utf-8")
    html = html.replace("__PRICER__", (ICI / "pricer.js").read_text(encoding="utf-8"))
    # quantiles.js lit Pricer au chargement : il vient APRES
    html = html.replace("__QUANTILES__", (ICI / "quantiles.js").read_text(encoding="utf-8"))
    html = html.replace("__BACKTEST__", (ICI / "backtest.js").read_text(encoding="utf-8"))
    html = html.replace("__MOTEURBT__", (ICI / "onglet_backtest.js").read_text(encoding="utf-8"))
    # la vue se branche sur le classeur : elle doit etre chargee en DERNIER
    html = html.replace("__VUEBT__", (ICI / "vue_backtest.js").read_text(encoding="utf-8"))
    html = html.replace("__FUTURES__", futures_json())
    html = html.replace("__DATA__", jour_json(jour0))
    html = html.replace("__JOURS__", json.dumps(index))
    # le gabarit vient de « CGMY a la main » : ici le modele est un CHOIX, pas
    # une identite. Le nom doit suivre, sinon la page ment sur ce qu'elle fait.
    html = html.replace("CGMY à la main", "Calibration à la main")

    cible = RACINE / "index.html"
    tmp = cible.with_suffix(".tmp")
    io.open(tmp, "w", encoding="utf-8").write(html)
    os.replace(tmp, cible)
    print("index.html : {} jours ({} integre), {:.1f} Mo".format(
        len(index), jour0, cible.stat().st_size / 1e6))


if __name__ == "__main__":
    main()
