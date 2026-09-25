# Calibration à la main — socle autonome

Dossier conçu pour être copié tel quel sur un poste isolé. **Aucune installation** :
ni Python, ni serveur, ni extension. Un navigateur Chrome ou Edge suffit.

## État

| | |
|---|---|
| `sources/pricer.js` | **prêt et validé** — le moteur complet |
| `sources/test_pricer.js` | **prêt** — la validation, rejouable par `node sources/test_pricer.js` |
| `jours/` | **prêt** — 253 instantanés quotidiens, 14/01/2024 → 24/09/2026 |
| `futures.js` | **prêt** — historique complet des futures, 286 échéances |
| `sources/telecharger_futures.py` | **prêt** — rafraîchit `futures.js` (poste connecté uniquement) |
| `index.html` | **à finir** — voir « Ce qui reste » |

## Le moteur

Un modèle se **compose** de blocs. L'exposant de la fonction caractéristique est
leur somme, et un bloc absent ne laisse **aucun paramètre** derrière lui :
Black-Scholes n'est pas un Heston à σ nul, c'est un exposant qui ne contient pas κ.

```
log Φ(u) = A_variance(u,T) + A_sauts(u,T)
```

**Bloc variance** — `p.var`

| type | paramètres | |
|---|---|---|
| `const` | `sig` | Black-Scholes |
| `heston` | `kappa, theta, sigma, rho, v0` | un facteur CIR |
| `heston2` | `f: [{…}, {…}]` | deux facteurs indépendants, exposants additifs |
| `hestonJ` | `… , lamV, muV` | + sauts de variance **exponentiels** |

**Bloc sauts du spot** — une clé au choix, ou aucune

| clé | paramètres |
|---|---|
| `sauts` | `[{lam, muJ, sigJ}, …]` — 1 composante = Merton, 2 = bi-lognormal |
| `kou` | `{lam, p, eta1, eta2}` — double exponentielle, **η₁ > 1 obligatoire** |
| `cgmy` | `{C, G, M, Y}` — **M > 1 obligatoire** |
| `vg` | `{sigma, nu, theta}` — Variance Gamma |
| `nig` | `{alpha, beta, delta}` — **α > \|β+1\| obligatoire** |

Les contraintes marquées « obligatoire » ne sont pas décoratives : en dessous,
`E[e^J]` diverge, le compensateur n'existe pas et le modèle n'est plus une
martingale. Les bornes doivent les porter.

**4 × 7 = 28 compositions**, toutes pricées par COS et par PROJ.

## Validation

`node sources/test_pricer.js` — trois contrôles :

```
1. MARTINGALE     call(K=F) − put(K=F)         0.0000 bps sur les 28
2. BLACK-SCHOLES  σ = 20/45/90 % → IV identique  < 0.005 bps
3. COS contre PROJ  écart max sur un smile     0.0000 bps
```

Le premier attrape un compensateur faux — la faute la plus facile en ajoutant une
famille. Le deuxième valide toute la chaîne, du CF à l'inversion Black-76. Le
troisième compare deux méthodes numériques indépendantes.

## L'instantané du jour, sans Python

La page (une fois finie) appellera Deribit directement :

```
public/get_book_summary_by_currency?currency=BTC&kind=option
public/get_index_price?index_name=btc_usd
```

Elle reconstruit `T`, `F` par échéance, `r = log(F/X)/T`, et inverse les IV
bid/ask/mid en Black-76 avec le pricer embarqué — le format exact de
`jours/<date>.js`.

**Pourquoi ça marche depuis un simple fichier.** Deribit renvoie
`Access-Control-Allow-Origin` en réfléchissant l'origine. Une page ouverte en
`file://` envoie `Origin: null`, et Deribit répond `null` — testé. Pas de proxy,
pas de serveur. Le seul prérequis est un accès réseau à `deribit.com` ; sans lui,
le bouton échoue proprement et les 253 jours embarqués restent disponibles.

## Les futures

`futures.js` porte l'historique de **286 échéances** : 118 futures réellement
cotés chez Deribit, repris de leur émission à leur expiration, et 168 forwards
implicites relevés dans les instantanés — Deribit ne cote pas de future à
échéance quotidienne, ces expirations-là n'ont pas d'instrument.

`sources/telecharger_futures.py` le rafraîchit, en incrémental : une échéance
expirée ne bouge plus et n'est jamais retéléchargée. À lancer depuis un poste
connecté, puis recopier `futures.js`.

## Ce qui reste

`index.html` — l'interface. Le moteur est prêt et vérifié ; il manque la page qui
l'expose. Elle reprendra celle de `cgmy_a_la_main` (onglets, empilement des jours,
paramètres propres à chaque jour, filtres à seuil de maturité, métriques,
optimiseur local avec sa garde), **moins** l'évolution différentielle, les
surrogates v5/rawy et l'écouteur Python, **plus** :

- deux sélecteurs de composition — bloc variance, loi des sauts ;
- le bouton « Récupérer le snapshot du jour » ;
- un séparateur d'année dans la liste des jours (`──── 2024 ────`) ;
- l'écriture directe des résultats dans `resultats/`, un fichier par
  **modèle × jour × jeu d'échéances × perte × bornes**.

Je ne l'ai pas livrée à moitié faite : une page de cette taille non testée dans
un navigateur casse, et cette semaine a montré ce que coûte une régression
silencieuse dans le pricer.
