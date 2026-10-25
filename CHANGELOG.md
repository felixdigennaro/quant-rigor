# Changelog

## v1.6.2, check per qconfluence e qregime, e due difetti che hanno trovato (2026-10-25)

Fino a qui nessun check toccava `qconfluence` e `qregime`. La suite ne ha 45 nuovi, tutti
contro riferimenti indipendenti: forme chiuse calcolate nel test, effetti piantati nei dati, una
scansione delle finestre rifatta con la sola libreria standard, il controllo di potenza del
modulo. I due moduli importano numpy, quindi senza numpy le loro sezioni compaiono come `SKIP`.

### Correzione: il cancello 2 di `confluence_gate` dipendeva dalla media della strategia base
Confrontava il t della pendenza, che non vede la media della strategia, con il nullo di
`search_null`, che misura il t della media di ogni regola contro zero. Su 40 serie per riga, un
effetto vero passava il cancello nello 0% dei casi con media base 1.0, e il rumore puro nel
100% con media -0.3. Ora il cancello centra il P&L prima di chiamare `search_null`: l'effetto
vero passa nel 100% dei casi con media 0, 0.375 e 1.0, il rumore nello 0% con 0, -0.3 e +0.3.
`search_null` non cambia; la sua docstring e §5b dicono di centrare quando lo si usa da solo.

### Correzione: `_spearman` sui pareggi
Usava ranghi ordinali, e un conteggio di confluenze è quasi tutto pareggi: su un conteggio 0-3
dava 0.2479 contro 0.2623 del manuale. Ora usa i ranghi medi. Il campo è solo informativo,
nessun cancello lo legge.

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
