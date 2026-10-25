# Changelog

## v2.0.0, i cancelli 2 e 5 rifatti (2026-10-25)

**Cambi che rompono.**
- Il cancello 2 di `confluence_gate` non gira più con le sole maschere. Serve
  `build_score(masks, pnl)`, la funzione che ha costruito il punteggio, oppure
  `preregistered=True` se il punteggio era fissato prima di guardare questo P&L. Senza nessuno
  dei due il cancello 2 non è eseguito e il verdetto è NON DECIDIBILE; i due insieme sollevano
  `ValueError`.
- `min_n` sparisce dalla firma di `confluence_gate`: serviva solo a `search_null`, che resta
  com'era.
- `res["search_null"]` cambia forma. Non ha più `ok`, `max` e `n_rules`; `q95` e `median` sono
  quantili di |t| del nullo. Si aggiungono `method`, `two_sided`, `block`, `n_rep`,
  `real_abs_t`, `p`, `null_abs_t` (un array numpy) e, con `build_score`, `n_drawn` e
  `n_rebuilt_constant`.
- Il cancello 5 esamina il bucket in basso quando la relazione scende e boccia un bucket che
  perde soldi. `res["concentration"]` aggiunge `bucket`, `bucket_n` e `bucket_mean`.

### Il cancello 2
- **Rifà la ricerca.** Mescola il P&L a blocchi, chiama `build_score` su copie del P&L mescolato
  e delle maschere, e ricalcola lo stesso t. Il punteggio vero deve essere quello che
  `build_score(masks, pnl)` ricostruisce. Una funzione che restituisce lo stesso punteggio su
  ogni P&L mescolato non dipende dal P&L, quindi non rifà la ricerca: il cancello non gira e la
  nota rimanda a `preregistered=True`.
- **Un t robusto.** La pendenza è studentizzata con l'errore standard HC3 (DiCiccio e Romano
  2017). Con il t OLS, rumore la cui varianza cresce col punteggio passava 78 volte su 400; con
  HC3 23.
- **Un p-value esatto, a due code.** Passa se (1 + #{|t nullo| ≥ |t reale|}) / (1 + n) ≤ 0.05.
  Il confronto con il 95° percentile interpolato era troppo largo con poche estrazioni (il 9.2%
  del rumore con 20, misura della revisione indipendente). Sotto 19 estrazioni il p-value non
  può arrivare a 0.05 e il cancello non gira. A due code, una confluenza che peggiora il P&L
  può passare, come già nei cancelli 1, 3 e 4.
- **Solo le estrazioni in cui la ricerca trova qualcosa.** Su un P&L mescolato una ricerca a
  soglia spesso non trova niente. Quelle estrazioni non entrano nel nullo: il cancello ne
  estrae finché `n_rep` trovano qualcosa, entro 10·`n_rep`, e ne vuole almeno 19. Contarle come
  |t| = 0 faceva passare il rumore proprio nei dataset in cui la ricerca vera aveva trovato
  qualcosa, cioè gli unici su cui qualcuno esegue il cancello. Con lo script della revisione
  indipendente e la ricerca a soglia della suite, il rumore passava il 9.6%, il 38.6% e il
  100% delle esecuzioni a soglia 2.5, 3 e 3.5; sul codice finale il 3.9%, il 4.0% e il 5.2%.
- **Blocchi.** Mescola blocchi di 25 periodi (`block`); `block=1` mescola i singoli periodi. Se
  più del 5% delle estrazioni pareggia il t reale il cancello non può passare qualunque siano
  i dati, quindi non gira. Succede con 2 o 3 blocchi, dove la serie stessa è una permutazione
  su 2 o su 6, e con un punteggio periodico il cui periodo divide `block`: il giorno della
  settimana con blocchi di 25 giorni dava p 1.0 e NON USARE su un effetto con |t| 20.
- Valori non finiti, punteggio o P&L costanti, un t non definito: il cancello non gira, con
  una nota che dice quale dei casi è.

### Il cancello 5
Esaminava sempre il bucket col punteggio più alto. Con il cancello 2 a due code una relazione
che scende ci arrivava, e il cancello 5 esaminava il bucket peggiore: col suo totale negativo
`top3_share` usciva negativo e il cancello passava sempre. Sugli stessi dati, con i guadagni
del setup migliore concentrati in 3 periodi, il verdetto era UTILIZZABILE 0 volte su 50 col
punteggio che sale e 49 su 50 col punteggio rovesciato; ora 0 su 50 in entrambi i versi.
Esamina il bucket in alto se la pendenza sale e quello in basso se scende. Quel bucket deve guadagnare (media > 0), restare
dello stesso segno a ogni cancellazione e non dovere metà del totale a tre periodi. Così
boccia anche una strategia che perde in ogni bucket, che prima poteva passare.

### Misurato
Passaggi del cancello 2 su 200 serie per riga, con 200 permutazioni e il `confluence_gate` vero
delle due versioni. Il rumore dovrebbe passare 10 volte su 200.

| forma | caso | v1.6.2 | v2.0.0 |
|---|---|---|---|
| punteggio 0-3 fissato prima | rumore | 3 | 9 |
| punteggio 0-3 fissato prima | effetto +0.06 / -0.06 (la confluenza peggiora) | 161 / 0 | 165 / 171 |
| 3 regole migliori su 200, cercate | rumore, media base 0 / 1.0 / -0.3 | 152 / 78 / 142 | 14 / 16 / 13 |
| 3 regole migliori su 200, cercate | effetto 0.15 sulle 3 regole | 195 | 140 |
| una regola che copre 25% / 50% / 75% | rumore | 8 / 14 / 19 | 6 / 6 / 6 |
| una regola che copre 50% | effetto 0.12 | 164 | 130 |
| P&L a due valori | rumore / effetto 0.25 | 0 / 0 | 6 / 176 |
| P&L e punteggio persistenti | rumore / effetto 0.1 | 26 / 192 | 7 / 150 |
| varianza che cresce col punteggio | rumore / effetto 0.1 | 11 / 146 | 10 / 120 |
| rendimenti futuri sovrapposti, 20 periodi | rumore, block 25 / 100 | 39 | 22 / 14 |
| rendimenti futuri sovrapposti, 60 periodi | rumore, block 25 / 100 | 60 | 34 / 12 |
| andamento comune nel tempo | rumore | 185 | 188 |

Il punteggio fissato prima è una riga sola perché il t della pendenza non vede la media della
strategia base. Dove la v2.0.0 passa meno effetti veri, la v1.6.2 li passava insieme al rumore
(fino al 76% sul punteggio cercato, il 13% sulle serie persistenti), oppure la v2.0.0 paga le
due code e l'errore standard robusto. Sul punteggio cercato il test è esatto: la revisione
indipendente misura 5.15% ± 0.55 su 6000 serie, con p-value uniformi.

### Limiti che restano
- **Una memoria del P&L più lunga del blocco.** Con rendimenti futuri sovrapposti su 60 periodi
  il rumore passa 34 volte su 200 con blocchi di 25 e 12 con blocchi di 100. Il blocco deve
  coprire la memoria del P&L, almeno l'orizzonte dei rendimenti sovrapposti.
- **Un andamento nel tempo comune a punteggio e P&L** passa in 188 casi su 200: nessun
  rimescolamento lo distingue da un effetto. Lo devono fermare i cancelli 3 e 4.
- **Una ricerca rifatta a metà.** `build_score` deve contenere ogni passo che ha guardato il
  P&L. Con una lista corta di 10 regole su 200 fissata sul P&L vero e solo la scelta delle 3
  rifatta, il rumore passa il 69% delle volte, contro il 4.75% rifacendo entrambi i passi
  (misura della revisione indipendente). Il codice non può accorgersene.
- **Un punteggio a code pesanti con varianza che cresce col punteggio.** HC3 non basta: con un
  punteggio lognormale il rumore passa 82 volte su 1000. Passando i ranghi del punteggio,
  44 su 1000.
- **Il cancello 5 guarda un bucket, non la pendenza.** Una pendenza dovuta a pochi periodi
  nell'altro bucket gli sfugge: tre crolli di -60 nel bucket che il dimensionamento pesa di
  meno portano il verdetto a UTILIZZABILE 33 volte su 100, contro 5 senza crolli (misura della
  revisione indipendente). Il rimedio proposto, il jackknife della pendenza invece che di un
  bucket, non è stato fatto.
- **Costo.** `build_score` gira fino a 10·`n_rep` volte, e `null_abs_t` e la matrice delle
  permutazioni crescono con `n_rep` × n.

### Check
La suite passa da 424 a 471 check. 48 mutazioni del codice nuovo rendono rossa la suite,
nessuna sopravvive.

## v1.6.2, check per qconfluence e qregime, e due difetti che hanno trovato (2026-10-25)

Fino a qui nessun check toccava `qconfluence` e `qregime`. La suite ne ha 45 nuovi. Quasi tutti
confrontano con un riferimento indipendente: forme chiuse calcolate nel test, effetti piantati
nei dati, una scansione delle finestre rifatta con la sola libreria standard. Quattro fanno
girare il controllo di potenza del modulo (`validate_power`) e misurano quanto spesso scatta sul
rumore puro, su 10 serie e non su una sola, perché su una serie scatta per caso circa una volta
su dieci. I due moduli importano numpy, quindi senza numpy le loro sezioni compaiono come `SKIP`.
Ogni check è stato provato rompendo il codice che copre: 36 mutazioni, tutte prese, e ognuno dei
45 check diventa rosso sotto almeno una. La prima stesura ne aveva tre che non potevano scattare
(un detector simmetrico, un conteggio uniforme su cui anche i quantili funzionano, un `KeyError`
che faceva crollare la suite prima del riepilogo); sono stati riscritti.

### Correzione: il cancello 2 di `confluence_gate` dipendeva dalla media della strategia base
Confrontava il t della pendenza, che non vede la media della strategia, con il nullo di
`search_null`, che misura il t della media di ogni regola contro zero. Su 40 serie per riga, un
effetto vero passava il cancello nello 0% dei casi con media base 1.0, e il rumore puro nel
100% con media -0.3. Ora il cancello centra il P&L prima di chiamare `search_null`. Con un
punteggio fissato prima della ricerca, l'effetto vero passa nel 100% dei casi con media 0, 0.375
e 1.0, il rumore nello 0% con 0, -0.3 e +0.3. `search_null` non cambia; la sua docstring e §5b
dicono di centrare quando lo si usa da solo.

### Correzione: `_spearman` sui pareggi
Usava ranghi ordinali, e un conteggio di confluenze è quasi tutto pareggi. Sui 400 valori del
check (conteggio 0-3 sbilanciato, seme 11) dava 0.1685 contro 0.1927 dei ranghi medi. Ora usa i
ranghi medi. Il campo è solo informativo, nessun cancello lo legge.

### Limiti noti del cancello 2, misurati e non riparati
Centrare toglie la dipendenza dalla media della strategia base, non gli altri difetti di questo
nullo. Tutti vengono dalla stessa premessa: il nullo misura il miglior t della media di una
regola su un sottocampione, il cancello confronta il t della pendenza del punteggio sul campione
intero. Misure del 2026-10-25, 40 serie per caso salvo dove indicato:
- **Punteggio costruito dalla ricerca.** Senza effetto, punteggio = quante delle 3 regole migliori
  su 200 (scelte sulla parte di addestramento) sono attive: il cancello 2 passa 31 volte su 40.
  Il nullo non ripete la scelta. Lo ferma il cancello 4, e senza replica il verdetto resta NON
  DECIDIBILE.
- **Una regola larga.** La permutazione estrae la parte di addestramento senza reinserimento
  dalla serie centrata, quindi il t del nullo si restringe di circa √(1−f). Senza effetto, con
  una sola regola che copre il 25%, 50% o 75% dei periodi, il cancello 2 passa 2, 2 e 5 volte
  su 40.
- **P&L a due valori** (take profit e stop fissi) con regole di 30-36 periodi. Una regola
  permutata può avere solo vincite, il suo errore standard tocca il pavimento di 1e-12 e il suo
  t arriva a circa 1e12. Il 95° percentile del nullo ha mediana 9e11 e il cancello 2 passa 0
  volte su 20 anche con un effetto piantato: potenza nulla.
- **Andamenti nel tempo.** Un punteggio che sale col tempo e un P&L che deriva di 0.3 lungo il
  campione, senza legame diretto fra i due: il cancello 2 passa 35 volte su 40. La permutazione
  presume periodi scambiabili.

Il rimedio di fondo è un nullo che permuta e ricalcola lo stesso t della pendenza, ripetendo la
scelta del punteggio quando gli viene data. Cambia cosa significa il cancello, quindi è una
decisione, non una correzione.

### Non toccato
Il cancello 2 resta a una coda (`slope_t > q95`) mentre i cancelli 1, 3 e 4 accettano anche una
pendenza negativa. Su 40 serie con pendenza -0.25, e media base 0, 1.0 o -0.3, il cancello 1
passa 40 volte su 40 e il cancello 2 nessuna: una confluenza che peggiora il P&L non arriva mai
a «utilizzabile». Se debba poterci arrivare è una scelta, non un difetto.

## v1.6.1, la suite sa fallire (2026-10-25)

`tests/test_all.py` stampava `SUITE RED` e usciva comunque con 0, quindi un check di CI
costruito sopra non poteva mai diventare rosso. La copia nel repo `quant-rigor` (1.3.0) aveva
ancora `sys.exit(1 if FAIL else 0)`; la copia installata come skill l'aveva perso.

- **L'uscita è tornata**, e la suite controlla in testa che resti l'ultima istruzione del
  file. Se manca, se diventa `sys.exit(0)` o se un check le viene scritto sotto (dove non
  girerebbe mai), la suite esce con 1 prima di eseguire un solo check.
- **Senza numpy la sezione ndarray non fa più crollare la suite.** I moduli che la suite
  copre usano solo la libreria standard, quindi senza numpy la sezione compare come `SKIP` e
  il riepilogo la conta (`sections skipped`). Prima l'`import numpy` falliva dopo tutti gli
  altri check e il riepilogo non veniva mai stampato.
- **README e SKILL.md non scrivono più il numero di asserzioni a mano.** Dicevano 255 e 303;
  il 2026-10-25 la suite ne contava 379 con numpy e 374 senza. SKILL.md dice anche che
  `qconfluence` e `qregime` non hanno ancora check.
- **Le due copie sono di nuovo uguali.** Il repo riceve da quella installata le versioni
  v1.4-v1.6 e i test ndarray del 2026-08-10, che non avevano una voce qui.

## v1.6 — le affermazioni si generano, non si scrivono (2026-08-01)

Nata dal fallimento documentato delle regole in prosa. Cinque frasi consegnate contraddicevano
una tabella stampata **due righe sopra**; le quattro regole scritte per impedirlo (R10, R11,
G8, G9) sono state violate tutte entro 48 ore. Le guardie eseguibili, nello stesso periodo,
hanno preso ogni bug.

### Nuovo — `scripts/qclaim.py` (formule 65-67) e `SKILL.md §0b G10`
- `cmp_` — la frase di confronto generata dai valori: non può invertirsi per distrazione
- `rel` — verifica di soglia che scrive sé stessa, con `ClaimError` opzionale
- `verdict` — deriva il verdetto dai criteri **dichiarati prima**, invece di sceglierlo dopo
- `audit` — marca ogni numero della bozza non tracciabile a una quantità calcolata
- `pct_change` — variazione col segno giusto (niente "aumenta del 25%" su un calo)

### Correzione — `qportfolio.effective_breadth`
Avvertenza esplicita: `N/(1+(N-1)ρ)` è la dimensione efficace per la **varianza di una media**;
trasportarla alla soglia sul **massimo** in multiple testing è un'euristica, non un risultato.
Per la penalità da ricerca si usa RAS, che la misura invece di dedurla da un ρ ipotizzato.
*(Correzione arrivata da una sessione parallela che ha ricalcolato invece di fidarsi.)*

### Test
343 → **374 asserzioni, tutte verdi**. Ogni nuovo test riproduce una frase realmente
consegnata e sbagliata.


## v1.5 — le guardie contro il fallimento silenzioso (2026-07-30)

Nata da una sessione in cui **nove errori miei** sono arrivati fino all'output, e nessuno di
essi ha prodotto un errore a runtime. Tutti erano della stessa famiglia: un controllo che non
puo' scattare, una regola che non puo' sparare, un numero che non ho misurato io.

**Chiude il buco complementare a v1.4:** v1.4 sapeva misurare il costo; v1.5 sa accorgersi
quando la misura non c'e' stata.

### Nuovo — sezione `SKILL.md §0b`, nove guardie in testa all'ordine delle operazioni
G1 soglia raggiungibile · G2 tasso di eventi · G3 scala log · G4 due guardie di plausibilita'
con applicabilita' · G5 il numero in prestito · G6 breadth con rho negativa · G7 il segno prima
della diversificazione · G8 hai scaricato quello che credi · G9 la prosa contro la tabella.

### Nuovo — `qstats.py` (formule 60-62)
- `reachable_threshold` — |z|max = `sqrt(W-1)` (pop) o `(W-1)/sqrt(W)` (camp.). Con W=4 una
  soglia |z|>2 e' irraggiungibile: la strategia stava in posizione lo 0.0% del tempo **senza
  errori**. Propone la finestra minima.
- `event_rate_check` — il tasso osservato contro quello inteso. Ha preso due regole
  "una al giorno" che ne facevano 7.5 (5 492 cambi/anno invece di 730).
- `log_to_simple` — conversione obbligatoria prima di pubblicare: −0.652 log = **−48%**, non −65%.

### Nuovo — `qmicro.py` (formula 63)
- `plausibility_flag` **riscritta**: seconda guardia fisica (`stima/range_barra > 0.5` e'
  impossibile) piu' l'**applicabilita'** della prima. Su BTC anche lo spread vero di 1bp
  implica 580 tick, quindi la guardia a tick produceva un **falso positivo**: ora si dichiara
  non applicabile e, se nessuna guardia si applica, la funzione **non emette giudizio**.
  Corretto anche un numero sbagliato nel docstring precedente (53 tick, non 178).
- `sign_flip_parameter` — bisezione sulla soglia a cui la conclusione cambia segno.
- `borrowed_parameter` — registra un parametro di terzi e **blocca la consegna del verdetto**
  se il ribaltamento cade nell'intervallo plausibile. Segnala le differenze di **categoria**,
  non solo di grandezza.

### Nuovo — `qportfolio.py` (formula 64)
- `effective_breadth` — **rifiuta** di chiamare breadth la formula quando `rho <= 0` (dava
  11.71 su 7) e riporta la firma della cancellazione: `rho <= 0` insieme a PC1 alta.
- `diversification_sign_check` — con `mu <= 0` diversificare **peggiora** lo Sharpe
  (misurato: 7 gambe a media −0.87 → aggregato −2.74).

### Test
307 → **343 asserzioni, tutte verdi**. Ogni nuovo test riproduce un errore realmente commesso.
Due dei nuovi test hanno **fallito alla prima esecuzione** e hanno scoperto la questione
dell'applicabilita' della guardia A — il caso in cui respingeva una stima corretta.


## v1.4 — microstructure, costs, correlated-search correction, mean reversion

Added after studying four public repositories (provenance in
`references/08-microstructure-costs.md`): machine-learning-for-trading,
Financial-Models-Numerical-Methods, awesome-quant, Computational-Finance-Course.

**Closes the biggest hole in the package: it could measure uncertainty but not cost.**

### New — `scripts/qmicro.py` (formulas 50-54)
- `roll_spread` — Roll (1984), with a *physical* guard: float rounding on a deterministic
  series produces alternating differences and thus a spuriously negative covariance, so a
  sign test alone is not enough.
- `corwin_schultz_spread` / `_median` — Corwin-Schultz (2012) high-low estimator; preferred
  on OHLC bars because it does not break on trending data.
- `bounce_check` — the test that settles short-horizon reversion: if the implied spread is of
  the order of one tick, the reversion **is** the cost of harvesting it.
- `cost_cliff`, `max_trade_frequency`, `required_edge_per_trade` — the frequency ceiling
  implied by a target net Sharpe. Run before building any intraday system.
- `decompose_slippage` — splits a gap between theoretical and realised P&L into entry slip,
  exit slip, commission and residual. Three components, three different fixes.

### New — `scripts/qpairs.py` (formulas 55-59)
- `ou_half_life`, `ou_fit` — half-life sets the trade count and therefore the cost hurdle;
  compute it before the backtest, not after.
- `hurst` — second moment about zero (the definition), not variance about the mean. Using the
  variance silently removes drift and returns H≈0 for a pure trend.
- `adf_stat` — returns the statistic and the critical values, deliberately **no p-value**:
  the null distribution is non-standard and a fabricated p-value is false precision.
- `engle_granger` — with the correct, more severe EG critical values (5%: −3.34 vs −2.86).
- `kalman_hedge_ratio` — time-varying beta, and a warning when beta drifts monotonically,
  which means the relationship is dying rather than oscillating.
- `pair_report` — one call, ending in estimated trades/year and annual cost.

### New — RAS in `scripts/qvalidate.py` (formula 38)
- `rademacher_complexity`, `massart_bound`, `effective_trials`, `ras_sharpe_adjustment`.
- Why: the Deflated Sharpe corrects for the *number* of trials assuming independence. A real
  parameter sweep produces highly *correlated* variants, so DSR over-penalises. Measured on
  40 synthetic variants over T=1000: correlated family → 1.1 effective trials; independent
  family → 10.4. Report `effective_trials`, not the raw configuration count.

### Protocol
- New §8b "COST IS A MEASUREMENT, NOT AN ASSUMPTION" — five steps before publishing any net
  figure; any unmeasured cost stays marked `⚠ assumed` and every dependent conclusion is
  stress-tested at 2×.

### Tests: 255 → 303 assertions, all green
Four new tests failed on first run. **The estimators were right; the tests were wrong** —
they asserted a property of a single realisation of a 5%-level test, which fails 5% of the
time by construction. Rewritten as calibration checks over many replications. Measured
afterwards: ADF false-positive rate 3.3% at the 5% threshold; half-life median 6.79 over 200
replications against a closed form of 6.58. **Rule now encoded: to test a statistical
procedure, test its operating characteristics, never a single draw.**

Versioning: MAJOR.MINOR.PATCH. MAJOR changes the protocol's steps or verdict semantics; MINOR adds
formulas, modules or functions; PATCH fixes correctness or presentation.

## 1.3.0 — 2026-07-26

### Added to `qdata.py`
- **Databento** support (paid, 16 years of CME/CBOT/NYMEX/COMEX plus ICE, equities, options):
  `db_cost()` estimates USD and record count BEFORE spending anything, and `db_download()`
  refuses to run without `confirm=True`. The API key is read from `DATABENTO_API_KEY` in the
  environment and is never accepted as an argument, printed, or written to disk.
- **Binance** free crypto: `binance_months()` / `binance_fetch()` pull the public monthly kline
  dumps — 1-minute bars from 2017-08, no API key, no rate limit. Verified reachable (HTTP 200,
  about 2 MB per month per symbol).
- **`CFD_TO_FUTURES` + `cfd_note()`**: a mapping from what CFD brokers call an instrument to the
  listed contract that actually has history. CFDs are OTC, so no vendor sells canonical CFD
  history — the honest method is to test on the underlying and add your broker's spread and
  financing as measured costs. 20 common CFD names mapped.

### Measured
- Databento historical pricing is $0.50/GB with $125 of signup credit, i.e. about 250 GB.
  NQ and GC both carry history from 2010-06-06 with OHLCV-1m available, so the datasets needed
  for an intraday study cost well under a dollar.

## 1.2.0 — 2026-07-26

### Added
- `references/07-montecarlo.md` — formulas 45-49.
- `qmontecarlo.py` — the null distributions a t-statistic assumes away: permutation and sign-flip
  tests (with the note that permuting order cannot test a mean), surrogate paths under four
  explicit no-edge processes (block/iid bootstrap, GBM, OU), White's Reality Check and Hansen's SPA
  for post-search comparison against a benchmark, and moment stress across six adversarial scenarios.
- `qdata.py` — the only networked module: public OHLCV history, your own CSV files, an idempotent
  forward collector, and `availability()` which MEASURES what a provider actually returns.
- 35 further assertions (255 total), including nulls validated against known cases: noise gives
  p > 0.10, a real edge p < 0.01, permutation ignores order on iid data but detects timing structure,
  Reality Check rejects a screen of eight noise strategies and finds a genuinely superior one.

### Documented (measured, not assumed)
- Free intraday history is the binding constraint on evaluating intraday strategies: ~20 years of
  daily bars and ~2 years hourly, but only ~60 calendar days of 5-minute (~50 sessions) and ~7 days
  of 1-minute. Older intraday windows return empty. A daily strategy is answerable from free data;
  an intraday one requires user-supplied history or forward collection. `SKILL.md` now checks data
  availability before agreeing to an analysis.

## 1.1.0 — 2026-07-26

### Added
- `qportfolio.py` — makes formulas 9-15 executable: covariance/correlation, portfolio volatility,
  marginal risk contributions, risk parity, Jacobi eigendecomposition, PC1 share, condition number,
  Marchenko-Pastur noise edge, shrinkage, capped long-only minimum variance.
- **`tail_dependence()`** — `λ = P(X ≤ q | Y ≤ q)`, the unbiased answer to "do these fail together".
- `qreport.py` — assembles an analysis into one self-contained HTML report: assumptions,
  measurability, full metric set with standard errors, paired difference, contraction check,
  overfitting diagnostics, sizing, five inlined SVG charts, verdict, and an explicit
  "what this evidence does not support" section. Objective register by construction.
- `examples/demo_portfolio.py` — third validation domain (four-sleeve book).
- 34 further assertions (220 total), including analytic portfolio identities (Euler risk
  decomposition, known eigenvalues, trace preservation, MP edge) and report-integrity checks.

### Fixed
- `solve_min_variance` violated the position cap: capping before renormalising let the
  renormalisation push weights back above it. Replaced with a water-filling projection.
- Correlation matrices could report a diagonal of 0 when an asset had zero variance in a subset;
  the diagonal is now 1 by convention.
- Documented a real statistical trap in `tail_conditional_corr`: Pearson correlation computed
  inside the worst decile of the portfolio mean is **biased downward** by the conditioning
  (Boyer-Gibson-Loretan), so a falling value is not evidence of improved diversification.
  `tail_dependence()` is now the reported measure; the biased one is kept only for comparison.

## 1.0.0 — 2026-07-26

First packaged release.

### Protocol
- Eight-step fixed order: restate, measurability gate, simulator validation, measure, uncertainty,
  contrarian pass, charts, verdict. Default action is "change nothing" absent decidable-positive evidence.
- Measurability gate promoted to step 2, before any backtest.

### Reference modules
- Formulas 1-30 across four modules (probability, linear algebra, stochastic calculus,
  pricing/risk), each entry with formula, practical use, **how it lies**, and a concrete action.
- Formulas 31-37 (`05-validation.md`): walk-forward, purged K-fold with embargo, Deflated Sharpe,
  PBO, HAC/Newey-West, MinTRL, multiple-testing/FDR.
- Formulas 38-44 (`06-sizing-execution.md`): Kelly and fractional Kelly, drawdown-constrained
  sizing, volatility targeting, market impact, capacity, turnover drag, latency and fill probability.
- `routing-table.md`: question-to-formula routing.
- `backtest-protocol.md`, `contrarian-protocol.md`, `charts.md`, `worked-example.md`.

### Code (stdlib only — no numpy, scipy, pandas or matplotlib)
- `qstats.py`: describe, paired difference, bootstrap (with block mode), jackknife, concentration,
  best-of-N headroom, required n, effective n, measurability gate, refinement factor, drawdown,
  empirical VaR/CVaR, Kelly ceiling.
- `qvalidate.py`: normal CDF/PPF, HAC standard errors, Deflated Sharpe, expected max Sharpe, PBO
  via CSCV, walk-forward and purged K-fold splits, MinTRL, Benjamini-Hochberg, Bonferroni.
- `qsize.py`: Kelly (continuous and discrete), growth rate, drawdown-constrained fraction,
  EWMA volatility, volatility targeting, square-root impact, capacity, turnover drag, cost ratio,
  effective rule distance under polling, adverse-fill adjustment.
- `qcharts.py`: equity, underwater, per-period contribution, parameter surface with noise band,
  return distribution with normal overlay — all as self-contained SVG. Plus ASCII sparkline.
- `tests/test_all.py`: 163 assertions against independent references (published t-tables, known
  normal quantiles, closed forms, analytic identities) and degenerate inputs.
- `examples/demo_strategy.py`: runnable end-to-end demo on a deterministic synthetic market.

### Fixed during development
- t-distribution tail was truncated, giving a 1% absolute p-value error at df=1 (Cauchy tails).
  Replaced with closed forms for df=1,2 and an analytic power-law tail remainder for df>=3;
  max error across 16 tabulated critical values now 3.3e-05.
- `dist_chart` ran off-canvas on a constant series (degenerate bin width).
- Zero-variance inputs produced ZeroDivisionError, 1e152 t-stats and a "strong evidence" verdict
  on NaN. All three now return explicit undefined states.
- `effective_n` could report more independent observations than the sample size under negative
  autocorrelation; now capped at n.
