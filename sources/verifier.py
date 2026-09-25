# -*- coding: utf-8 -*-
u"""Controle une page construite, avant de la croire bonne.

POURQUOI CET OUTIL. Les trois regressions de la journee avaient le meme
profil : une suppression trop fine qui coupe au milieu d'une balise, laisse un
fragment, et fait apparaitre le symptome AILLEURS -- un bouton qui disparait,
un panneau vide, un `id="r-de">` affiche en clair. `node --check` ne les voit
pas : le JavaScript reste valide, c'est le HTML qui ne l'est plus. Le navigateur
non plus ne dit rien : il recolle les morceaux en silence.

CE QU'ON CONTROLE
  1. l'equilibre des balises, par une vraie analyse syntaxique avec pile ;
  2. les fragments orphelins -- un attribut qui traine hors de toute balise ;
  3. la syntaxe de chaque bloc <script>, par node --check ;
  4. les marqueurs de gabarit non substitues (__PRICER__ et compagnie) ;
  5. la presence des identifiants dont le code depend.

Usage : python sources/verifier.py [chemin/vers/index.html ...]
        sans argument : ../index.html
"""
import io
import re
import subprocess
import sys
import tempfile
from html.parser import HTMLParser
from pathlib import Path

VIDES = {"area", "base", "br", "col", "embed", "hr", "img", "input", "link",
         "meta", "param", "source", "track", "wbr"}
# identifiants dont le code de la page se sert : leur absence est silencieuse
REQUIS = ["jour", "preset", "opt", "opt-prog", "status", "opt-meth", "feller",
          "dir-info", "reg-nom", "reg-save", "out", "meth", "cosN", "cosL",
          "mesure", "poids", "b-large", "b-prod", "excl", "dfotm",
          "bloc-var", "bloc-saut", "nempile", "emp-go", "parjour"]


class Pile(HTMLParser):
    u"""Verifie l'imbrication ; signale la premiere incoherence rencontree."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.pile = []
        self.soucis = []

    def handle_starttag(self, tag, attrs):
        if tag not in VIDES:
            self.pile.append((tag, self.getpos()[0]))

    def handle_startendtag(self, tag, attrs):
        pass

    def handle_endtag(self, tag):
        if tag in VIDES:
            return
        if not self.pile:
            self.soucis.append("ligne {} : </{}> sans ouverture".format(
                self.getpos()[0], tag))
            return
        if self.pile[-1][0] != tag:
            # on cherche plus bas : une balise a-t-elle ete laissee ouverte ?
            for k in range(len(self.pile) - 2, -1, -1):
                if self.pile[k][0] == tag:
                    manquants = [t for t, _ in self.pile[k + 1:]]
                    self.soucis.append(
                        "ligne {} : </{}> ferme, mais {} reste(nt) ouverte(s) (ouvert ligne {})"
                        .format(self.getpos()[0], tag, ", ".join("<%s>" % m for m in manquants),
                                self.pile[k][1]))
                    del self.pile[k:]
                    return
            self.soucis.append("ligne {} : </{}> inattendu, on attendait </{}>".format(
                self.getpos()[0], tag, self.pile[-1][0]))
            return
        self.pile.pop()


def scripts(html):
    for m in re.finditer(r"<script(?![^>]*\bsrc=)[^>]*>(.*?)</script>", html, re.S):
        yield m.group(1)


def controler(p):
    html = io.open(p, encoding="utf-8").read()
    print("\n=== {} · {:.1f} Mo ===".format(p.name, p.stat().st_size / 1e6))
    ko = 0

    # --- 1 et 2 : le balisage, scripts retires (ils contiennent des '<') ---
    sans = re.sub(r"<script.*?</script>", "", html, flags=re.S)
    sans = re.sub(r"<style.*?</style>", "", sans, flags=re.S)
    pile = Pile()
    pile.feed(sans)
    if pile.soucis:
        ko += len(pile.soucis)
        print("  BALISAGE :")
        for x in pile.soucis[:6]:
            print("    " + x)
    else:
        print("  balisage : equilibre")
    ouverts = [t for t, _ in pile.pile if t not in ("html", "body", "head")]
    if ouverts:
        ko += 1
        print("    restent ouvertes : " + ", ".join("<%s>" % t for t in ouverts[:8]))

    # fragment orphelin : un attribut hors de toute balise, comme id="r-de">
    orph = re.findall(r'^\s*(?:id|class|style|type|value)="[^"]*">', sans, re.M)
    if orph:
        ko += len(orph)
        print("  FRAGMENTS ORPHELINS : " + " · ".join(o.strip() for o in orph[:4]))
    else:
        print("  fragments orphelins : aucun")

    # --- 3 : la syntaxe des scripts ---
    mauvais = []
    with tempfile.TemporaryDirectory() as d:
        for i, code in enumerate(scripts(html)):
            f = Path(d) / "b{}.js".format(i)
            f.write_text(code, encoding="utf-8")
            r = subprocess.run(["node", "--check", str(f)],
                               capture_output=True, text=True, shell=True)
            if r.returncode:
                mauvais.append((i, (r.stderr or "").strip().split("\n")[:3]))
    if mauvais:
        ko += len(mauvais)
        print("  SCRIPTS :")
        for i, msg in mauvais:
            print("    bloc {} : {}".format(i, " / ".join(msg)))
    else:
        print("  scripts : {} bloc(s), syntaxe correcte".format(len(list(scripts(html)))))

    # --- 4 : marqueurs non substitues ---
    restes = sorted(set(re.findall(r"__[A-Z_]{3,}__", html)))
    if restes:
        ko += len(restes)
        print("  MARQUEURS NON REMPLACES : " + ", ".join(restes))
    else:
        print("  marqueurs : tous substitues")

    # --- 5 : identifiants attendus ---
    ids = set(re.findall(r'id="([^"]+)"', html))
    absents = [k for k in REQUIS if k not in ids]
    if absents:
        print("  identifiants absents : " + ", ".join(absents)
              + "   (normal si la page n'a pas cette fonction)")
    else:
        print("  identifiants : les {} attendus sont la".format(len(REQUIS)))

    print("  ==> {}".format("CONFORME" if ko == 0 else "{} PROBLEME(S)".format(ko)))
    return ko


def main():
    cibles = [Path(x) for x in sys.argv[1:]] or [Path(__file__).resolve().parent.parent / "index.html"]
    total = sum(controler(p) for p in cibles if p.exists())
    print("\n{}".format("Tout est conforme." if total == 0
                        else "{} probleme(s) au total.".format(total)))
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
