---
name: quant-rigor
description: Rigorous quantitative analysis of any trading or investment strategy — existing, proposed, or under evaluation. Applies 67 quant formulas (probability, linear algebra, stochastic calculus, pricing/risk, validation, sizing/execution, microstructure/costs, mean reversion, guards against silent failure) as working tools rather than decoration, each with its failure mode. Runs a measurability gate before any backtest, validates the simulator against reality, quantifies uncertainty (SE, t, bootstrap, jackknife, best-of-N headroom, Deflated Sharpe, RAS/Rademacher for correlated searches, PBO, HAC), estimates transaction costs from your own bars instead of assuming them, refuses to deliver a verdict that hinges on a number it has not measured itself, runs an adversarial contrarian pass, and produces the charts that expose what tables hide. Use whenever a strategy, edge, backtest result, parameter change, sizing decision, or performance claim is on the table.
---

# Quant Rigor — analyse before you believe

> A number without its standard error is an opinion. A backtest without an out-of-sample is a story.
> A strategy that "works" on 15 days works on 15 days.

This skill turns the 30 core quant formulas into an operating procedure. It is not a formula sheet
— the formulas are in `references/`. This file is the **workflow**: the order of operations that
stops a plausible-but-wrong result from reaching a decision that moves money.

---

## 0. WHEN THIS FIRES

Fires on any of these, automatically, without being asked:
- A strategy is proposed, changed, or parameter-tuned ("what if we add a trailing stop / move the TP / raise size")
- A backtest result is reported or quoted (mine or someone else's)
- A performance claim is made ("this returns X%", "the edge is real", "win rate is Y")
- A sizing / leverage / risk-budget decision is on the table
- An existing live strategy is reviewed
- Someone asks whether an edge is real

Does NOT fire on: pure execution plumbing (order routing, API wiring) with no claim about returns,
or accounting of realised P&L already reconciled to an exchange balance (that's bookkeeping, not inference).

---

## 0b. LE GUARDIE CONTRO IL FALLIMENTO SILENZIOSO — prima di tutto il resto

Il resto di questo protocollo difende da **numeri sbagliati**. Questa sezione difende dalla
classe che non produce numeri sbagliati affatto: **controlli che non possono scattare, regole
che non possono sparare, verdetti che poggiano su un numero preso in prestito.** Non danno
errore, non sembrano rotti, e assomigliano a risultati.

Ogni guardia qui esiste perche' l'errore corrispondente e' stato commesso davvero, e ogni
guardia ha un test che lo riproduce (`tests/test_all.py`, sezione v1.5).

### G1 · La soglia e' raggiungibile? — `qstats.reachable_threshold`
Uno z-score su finestra `W` ha un massimo **algebrico**: `sqrt(W-1)` con sd di popolazione,
`(W-1)/sqrt(W)` con sd campionaria. Con `W=4` il massimo e' **1.73**: una regola `|z|>2` non
puo' sparare **mai**. Risultato osservato: 0.0% del tempo in posizione, zero errori.
> **Un test che non puo' scattare non e' un test superato: e' un test rotto, e somiglia a un
> risultato negativo.** Chiamala prima di ogni backtest con una soglia su una statistica
> normalizzata.

### G2 · Il tasso di eventi corrisponde all'intenzione? — `qstats.event_rate_check`
Due regole "una operazione al giorno" ne facevano 7.5 al giorno perche' mancava il flag
"gia' fatto": 5 492 e 2 743 cambi/anno invece di ~730. Nessun errore, solo tutto il P&L a
valle sbagliato. Chiamala **subito dopo aver generato i segnali e prima di guardare il P&L**.
Zero eventi non e' un risultato negativo: e' una regola che non spara (→ G1).

### G3 · La scala del numero che pubblichi — `qstats.log_to_simple`
Drawdown calcolati su `cumsum` di log-rendimenti e stampati come percentuali: `−0.652` in log
e' **−48%** reale, non −65%. Diciassette punti di errore su ogni riga di ogni tabella.
> Qualunque metrica nata in spazio log si converte **prima** di pubblicarla, e la scala si
> dichiara nell'etichetta.

### G4 · Due guardie di plausibilita', e sapere quando una non si applica
`qmicro.plausibility_flag` ora ne ha due, che coprono casi **disgiunti**:
- **A · tick impliciti** — prende gli stimatori di microstruttura su barre grosse (ES daily:
  62bp = 53 tick, respinta).
- **B · quota del range di barra** — `stima / range_mediano > 0.5` e' **fisicamente
  impossibile**, perche' ogni scambio che attraversa lo spread crea un range almeno pari.
  Prende i casi in cui il tick e' molto piu' fine dello spread (BTC: 17.87bp = 2.9x il range).

E la lezione che i test hanno tirato fuori dopo: **la guardia A va dichiarata non applicabile
quando il tick e' troppo fine.** Su BTC anche lo spread VERO di 1bp implica 580 tick, quindi
applicare A produce un **falso positivo** — respinge una stima corretta. Se nessuna guardia e'
applicabile la funzione **non emette giudizio**, invece di emetterne uno sbagliato.

### G5 · Il numero in prestito — `qmicro.borrowed_parameter`, `qmicro.sign_flip_parameter`
**La lezione piu' costosa.** Un penalty di −9c preso dal README di un terzo decideva un
verdetto. Misurato nel proprio contesto: **−44.3c, 4.9x piu' grande**, e per un **errore di
categoria** — le sue posizioni erano chiudibili, quelle in esame erano binari tenuti a
scadenza, senza uscita. Il verdetto si e' ribaltato.
> **Un numero preso in prestito non e' un numero verificato.** Non e' assunto (quindi R1 non
> lo prende) e non e' misurato da te. Regola:
> 1. calcola la **soglia di ribaltamento** del parametro (`sign_flip_parameter`);
> 2. se il ribaltamento cade **dentro** l'intervallo plausibile, **non consegnare il verdetto**
>    — consegna la soglia, e vai a misurare;
> 3. prima di prendere in prestito, cerca la **differenza di categoria** fra i due contesti,
>    non solo quella di grandezza.

### G6 · Breadth con rho negativa — `qportfolio.effective_breadth`
`N/(1+(N-1)rho)` applicata a 7 strategie sullo stesso strumento con `rho = −0.067` ha dato
**11.71 su 7**. La formula misura **riduzione di varianza**, non "scommesse indipendenti con
rendimento atteso": due strategie che si **annullano** hanno `rho = −1` e breadth infinita
mentre pagano due commissioni per tenere zero. La funzione ora **rifiuta** di chiamarla breadth
quando `rho <= 0` e riporta la **firma della cancellazione**: `rho <= 0` **insieme a** PC1 alta.
In quel caso la domanda giusta non e' la breadth ma: quanta esposizione **lorda** sopravvive
nel netto? (I costi si sommano, i rendimenti si annullano.)

### G7 · Il segno prima della diversificazione — `qportfolio.diversification_sign_check`
Sette strategie con Sharpe medio **−0.87** aggregate danno **−2.74**, peggio della media delle
parti. Non e' un paradosso: diversificare riduce il **denominatore**, e con numeratore negativo
amplifica il segno. **Le tecniche di riduzione della varianza si applicano solo dopo aver
risolto il segno di mu.**

### G8 · Hai scaricato quello che credi? (vale per qualunque fonte esterna)
Una query ordinata per data decrescente ha restituito mercati **a 24 ore nel futuro**; i loro
book erano segnaposto seminati dal market maker, e quattro minuti di polling stavano per
diventare la conclusione "il libro e' statico". Prima di **qualunque** inferenza su dati
esterni: **timestamp contro adesso, volume/trade diversi da zero, e almeno un caso di
controllo noto.** Un guasto di raccolta non e' un verdetto sul dato (R8).

### G9 · La prosa non puo' contraddire la tabella
Due volte nella stessa sessione il commento ha affermato il contrario del numero stampato due
righe sopra ("11.71 resta sotto 8.37"; "salva il 60-70%" con 25% misurato).
> Ogni numero nella prosa **si genera dalla stessa variabile** che ha prodotto la tabella —
> non si ribatte a mano. E prima di consegnare, la prosa si rilegge **contro** le tabelle.

### G10 · La frase si genera dai numeri — `qclaim.py`
**La guardia nata dal fallimento delle altre.** Cinque volte fra luglio e agosto 2026 una frase
ha contraddetto una tabella stampata *due righe sopra*: `"11.71 resta sotto 8.37"` · `"salva il
60-70%"` (misurato 25%) · `"3.50 è la media"` (era 3.013) · `"aggiornata la casella"` (grep
diceva 0) · costo/range applicato a uno strumento binario. In ogni caso **il numero calcolato
era giusto e la frase sbagliata**, e la smentita era già visibile nell'output.

Le regole in prosa contro questo errore — R10, R11, G8, G9 — sono state violate **tutte** entro
48 ore dall'essere scritte. Le guardie eseguibili dello stesso periodo hanno preso ogni bug. La
differenza non è la qualità della regola: **una funzione scatta da sola, una frase no.**

- `cmp_(a, b, la, lb)` — genera la frase di confronto **dai valori**: non può invertirsi.
- `rel(v, op, soglia, label)` — verifica la soglia e la scrive; `raise_on_fail=True` per bloccare.
- `verdict(metriche, criteri)` — **deriva** il verdetto dai criteri dichiarati a inizio lavoro,
  invece di sceglierlo a fine analisi guardando i risultati.
- `audit(testo, valori)` — scorre la bozza e marca ogni numero **non tracciabile** a una
  quantità davvero calcolata: prende quelli ricordati o inventati.

Da usare al passo 3 del checkpoint R5 (*ri-verifica delle affermazioni*), che è il passo che
ricontrolla le frasi mentre il passo 2 ricontrolla il codice.


---

## 1. THE ORDER OF OPERATIONS (never reorder)

The sequence matters because each step can **kill the analysis**, and killing it early saves
the rest. Most wrong conclusions come from doing step 4 before steps 1-3.

```
0. GUARDIE       G1-G10 (§0b): la regola puo' sparare? il numero e' tuo? la frase e' GENERATA?
1. RESTATE       what exactly is claimed, and what would falsify it
1b. DECIDABILITY  is the horizon long enough for this Sharpe to reach t=2? (§1b — cheapest gate)
2. MEASURABILITY can the available data even measure this? (§2 — the step everyone skips)
3. VALIDATE      does the simulator reproduce reality where reality is known? (§3)
4. MEASURE       run it, with costs, on the full sample (§4)
5. UNCERTAINTY   band not point: SE, t, bootstrap, jackknife (§5)
6. CONTRARIAN    adversarial pass whose mandate is to destroy the finding (§6)
7. CHARTS        the four plots that expose what the table hides (§7)
8. VERDICT       decidable / not decidable, and what would make it decidable (§8)
```

Load `references/` files as needed — **start with `routing-table.md`** (question → formula):
- `routing-table.md` — which formula answers which question, by stage of work
- `01-probability.md` — formulas 1-8 (Bayes → LLN): inference, what a sample can support
- `02-linear-algebra.md` — 9-15 (matrices → Lagrange): portfolios, risk decomposition, optimisation
- `03-stochastic.md` — 16-22 (Taylor → SDE): price dynamics, why paths behave as they do
- `04-finance.md` — 23-30 (PV → VaR): discounting, options, Sharpe, CAPM, tail risk
- `05-validation.md` — 31-37: walk-forward, purged K-fold, Deflated Sharpe, PBO, HAC, MinTRL, FDR
- `06-sizing-execution.md` — 38-44: Kelly, drawdown-constrained sizing, vol targeting, impact,
  capacity, turnover drag, latency/fills
- `07-montecarlo.md` — 45-49: permutation & sign-flip, surrogate paths, White's Reality Check,
  Hansen's SPA, moment stress
- `08-microstructure-costs.md` — 50-59: spread estimation (Roll, Corwin-Schultz), the cost cliff,
  slippage decomposition, Rademacher/RAS, half-life, cointegration, Kalman hedge ratio
- `backtest-protocol.md` — the full checklist, bias by bias
- `contrarian-protocol.md` — the attack playbook (run it, don't skim it)
- `charts.md` — the required plots and what each one is for
- `worked-example.md` — **a real case end to end**, with the four mistakes that shaped this
  checklist. Read this first if you want the skill in 5 minutes.
- `scripts/qstats.py` — describe / paired / bootstrap / jackknife / headroom / required_n /
  effective_n / measurability / refinement_factor / drawdown / VaR-CVaR / Kelly.
  `python3 qstats.py --a <base…> --b <variant…> --variants N` prints the whole decision set.
- `scripts/qvalidate.py` — HAC/Newey-West, Deflated Sharpe, expected-max Sharpe, PBO (CSCV),
  walk-forward & purged K-fold splits, MinTRL, Benjamini-Hochberg FDR.
- `scripts/qsize.py` — Kelly (continuous/discrete), drawdown-constrained fraction, EWMA vol &
  vol targeting, square-root impact, capacity, turnover drag, cost ratio, effective rule distance
  under polling, adverse-fill adjustment.
- `scripts/qportfolio.py` — covariance/correlation, **tail dependence** (the unbiased tail measure),
  portfolio vol, marginal risk contributions, risk parity, Jacobi eigendecomposition, PC1 share,
  condition number, Marchenko-Pastur noise edge, shrinkage, capped min-variance.
- `scripts/qreport.py` — assembles the whole analysis into ONE self-contained HTML report
  (assumptions · measurability · metrics with SE · difference · overfitting · sizing · charts ·
  verdict · what the evidence does not support). Objective register, no persuasive language.
- `scripts/qmicro.py` — **what it actually costs**, estimated from your own bars: Roll and
  Corwin-Schultz spread estimators, `bounce_check` (is that short-horizon reversion just the
  bid-ask bounce?), `cost_cliff` / `max_trade_frequency` (the frequency ceiling implied by a
  target net Sharpe), `decompose_slippage` (entry slip vs exit slip vs commission).
- `scripts/qpairs.py` — mean reversion & pairs: `ou_half_life` (compute it BEFORE the
  backtest — it sets the trade count and therefore the cost hurdle), `hurst`, `adf_stat`,
  `engle_granger` (with the correct, more severe EG critical values), `kalman_hedge_ratio`,
  `pair_report`.
- `scripts/qconfluence.py` — **conditional edges and confluence sizing** (§5b): `dose_response`
  (monotone gradient across buckets, not a binary filter), `search_null` (what t your search
  produces on shuffled data — the bar your real t must clear), `confluence_gate` (the five gates,
  verdict derived not chosen; gate 2 replays your search through `build_score`, or takes
  `preregistered=True`), `sizing_from_slope` + `sizing_cost` (weights from the fitted slope,
  and what sizing costs you in Sharpe if the confluence is noise).
- `scripts/qstats.py` also carries the **decidability gate** (§1b): `years_to_decide` (how many
  years this Sharpe needs to reach t=2, with and without publication decay) and
  `literature_queue` (sort published strategies by decidability — the stated Sharpe orders the
  queue, it is never evidence).
- `scripts/qregime.py` — **"no edge works forever"** (§5c): `window_scan` (strong windows vs a
  BLOCK-shuffled null — picking the best window is itself a search), `regime_switch` (ex-ante
  detector against always-on AND a random switch at equal duty cycle), `regime_report`,
  `validate_power` (proves the test can fire before you trust a negative).
- `scripts/qcharts.py` — stdlib-only, emits standalone **SVG** (no matplotlib):
  `equity_chart`, `drawdown_chart`, `contribution_chart`, `surface_chart(…, headroom=…)`,
  `dist_chart`, plus `sparkline()` for terminal output.

The scripts ship with `tests/test_all.py`, which checks them against independent references
(published t-tables, known normal quantiles, closed forms) plus degenerate-input cases. It prints
the passed and failed counts and exits 1 if any check fails. `qconfluence` and `qregime` import
numpy, so their checks run only where numpy is installed and are listed as SKIP elsewhere.
Run it first.

---

## 1b. DECIDABILITY GATE — is this answerable at all, before you write anything?

**The cheapest gate in the file, and it runs before measurability.** From `t = Sharpe·√years`:

```
years_needed = (t_target / Sharpe)²        qstats.years_to_decide()
```

A Sharpe of 0.35 needs **31 years** to reach t=2. There is no cleverness that fixes this — not
a better simulator, not more features, not a finer bar. If the horizon you have is shorter than
`years_needed`, the backtest cannot decide the question and writing it produces a number with
no inferential content.

### Reading published strategies — the queue, not the answer

Measured 2026-08-04 across 60 published strategies with stated Sharpe ratios
(`paperswithbacktest/awesome-systematic-trading`, which ships both the numbers and runnable code):

| | |
|---|---|
| median stated Sharpe | **0.359** → 31 years to t=2 |
| already **negative** in the paper | 9 of 60 |
| decidable in 20 years at the stated value | 22 of 60 |
| decidable in 20 years **after publication decay** | **0 of 60** |

Decay is not a guess: McLean & Pontiff (2016) measure ~58% of anomaly returns disappearing
out-of-sample across 97 anomalies — roughly half statistical noise in the original sample, half
genuine arbitrage after publication. `PUBLICATION_DECAY = 0.42` encodes the surviving fraction,
and it is itself a **borrowed number** (G5): it comes from US equity anomalies, so restate or
re-measure it when the strategy is of another kind.

> **A published Sharpe orders your queue. It is never evidence.** It is the maximum of a search
> you did not control, in a context you did not verify, filtered by a journal that does not
> publish nulls. Use `qstats.literature_queue()` to sort candidates by decidability and work
> down that order; never adopt the number itself as a premise (§9, G5).

Calibration that puts the table in perspective: **ES buy-and-hold over 20 years reaches t = +2.46**
— the most documented return in finance, barely clearing the bar. Anything claiming to beat it
needs to clear a bar that the benchmark itself only just clears.

For the rigorous version of this gate — correcting for skew, kurtosis and number of trials —
use `qvalidate.min_track_record_length` and `qvalidate.deflated_sharpe`. `years_to_decide` is
the one-line version, meant to stop work before it starts.

---

## 2. MEASURABILITY GATE — can the data measure the claim at all?

**Run this BEFORE any backtest. It is the single highest-yield step and the most commonly skipped.**

An exit rule, a signal, or a stop is only measurable if the data's resolution is **much finer**
than the distance the rule cares about. Otherwise you are measuring your own blindness.

### The test
1. Compute the **typical bar range** for your data: `median(high - low)` over the relevant session.
2. Compute the rule's **characteristic distance** `D` (trail width, stop distance, target size).
3. Ratio `R = D / median_range`.
   - `R ≥ 3-5` → measurable
   - `1 ≤ R < 3` → **not measurable**; results are artefacts
   - `R < 1` → the rule lives entirely inside one bar; the backtest is fiction

### The convergence check (decisive)
Re-run the identical rule at **two or more data resolutions** (e.g. 5m and 1m):
- **Fixed-level rules** (limit TP, fixed stop) → the number **converges**: same result at both
  resolutions, because a fixed level is either touched in the bar or not. Resolution only
  reorders two levels inside one bar (1 bit of information, and only when both are touched).
- **Path-dependent rules** (trailing stop, running-max/min, breakeven-move, anything defined
  ON the path) → the number **does not converge**: it keeps changing as you refine, because
  measured path length grows like `√(1/Δt)` (Wiener process — see `03-stochastic.md` §21).
  Measure the refinement factor: `sum(fine bar ranges) / coarse bar range`. For a random walk
  this is `√(n)`. If you measure ≈`√n`, the path has no intrinsic scale and **no resolution
  you have is enough** — only rules with `D` far above the finest bar range are testable.

### Do NOT do this
**Aggregating fine bars into coarse bars is not a granularity test.** Coarsening *hides* adverse
excursions, so it monotonically improves everything path-dependent — and eventually starts hiding
the baseline's stop too. It proves nothing. The informative direction is **refinement**, never coarsening.

### Verdict language
If `R < 3` and refinement doesn't converge, the honest output is:
> "Not measurable with this data. To answer it I need <resolution>, collected forward
> with a pre-registered criterion. What IS measurable now: <fixed-level alternatives>."

Offer the non-path-dependent alternative that answers the same *business* question
(e.g. instead of a trailing stop: partial TP, time-based exit, daily loss limit — all fixed rules).

---

### Data first — check before promising an answer
Run `qdata.availability(symbol)` (or `python3 scripts/qdata.py --symbol X`) BEFORE agreeing to
evaluate a strategy. Measured, not assumed: ~20 years daily, ~2 years hourly, ~60 calendar days of
5-minute, ~7 days of 1-minute. Intraday retention is a hard limit — older intraday windows return
empty. Therefore: a daily strategy is answerable from free data; an **intraday** one is not, unless
the user supplies history (`qdata.from_csv`) or a forward collector has been running
(`qdata.append_csv`). Combine with `qstats.required_n()` and say the number of periods needed out
loud. Never substitute a short intraday sample for an answer.

## 3. VALIDATE THE SIMULATOR AGAINST REALITY

Never trust a simulator you haven't tied to a known outcome.

- If there are **live fills**: reproduce them. Entry levels must match **to the tick**; realised
  P&L must reconcile within costs. Report the comparison explicitly.
- If there are no live fills: reproduce a **published/known** result of the same rule, or
  hand-compute 2-3 periods and diff.
- Reconcile the **cash**: `start balance + Σ realised P&L + funding − fees = end balance`.
  If it doesn't tie to the cent, stop — something in the model is not what you think.
- Then, and only then, believe its output on the untested variants.

Also validate the **direction of conversion**: points → currency → % of capital. A 10× error here
is invisible in the ratios and fatal in the sizing.

---

## 4. MEASURE — with costs, on everything

- **Full sample, no cherry-picking of periods.** Report `n` (independent observations, not rows).
- **Estimate the cost from your own bars first** (`qmicro`): Roll and Corwin-Schultz give a
  spread from OHLC alone, and `bounce_check` tells you whether a measured short-horizon
  reversion is anything more than the spread you would pay to chase it. Then:
- **Costs from reality**, never assumed: fee per round-trip and slippage per exit type,
  measured from actual fills. Express them as a **fraction of the target**:
  `cost / target`. If a rule's target is 40pt and round-trip cost is 7pt, that's **18%** — a rule
  whose cost ratio is 3-5× the baseline's is a different economic proposition, not a tweak.
- **Adverse selection on costs**: exits that fire during momentum bursts (stops, trailing) slip
  worse than calm scheduled exits (time flat). Do not reuse one slippage number for both.
- **Unmodelled costs must be named**: funding/carry per hour held, borrow, taxes (never model
  taxes unless that is the task), overnight gaps. If a variant holds longer, its number is
  optimistic by an unquantified amount — say so.
- **Executability**: is the tested parameter the one that can actually run? A trail of 40pt
  polled every 15s behaves like ~50-55pt. A stop triggered on mark/oracle is not a stop on last.
  Test what you can execute, not what you can imagine.

---

## 5. UNCERTAINTY — the step that decides most questions

Use `scripts/qstats.py`. Formulas in `01-probability.md` (§3-8).

### Minimum reporting set
- `mean ± SE` where `SE = σ/√n`, and `t = mean/SE`
- For a **difference** between two variants, use the **paired** difference per period
  (`SE_diff = σ_diff/√n`) — it removes the common market factor and is far more powerful
- `n`, and the **effective** n (clustered/overlapping observations count for less)
- **Bootstrap** the metric (10k resamples) → 5/95 percentile band
- **Jackknife**: drop one period at a time; does the sign survive all `n` deletions?
- Contribution concentration: what fraction of the result comes from the top 3 periods?

### Thresholds, stated plainly
| t | reading |
|---|---|
| < 1.0 | noise; do not act |
| 1.0-1.6 | suggestive; a hypothesis, not a finding |
| 1.6-2.0 | weak evidence (p ≈ 0.05-0.11); act only with a structural argument |
| > 2.0 | evidence, still conditional on one regime |
| > 3.0 | strong, if the sample spans regimes |

### Formal overfitting diagnostics (use them, not just the headroom)
- **Deflated Sharpe** (`qvalidate.deflated_sharpe`) — corrects the Sharpe for the number of trials,
  non-normality and sample length. `DSR < 0.95` → not established.
- **PBO** (`qvalidate.pbo`) — probability that the *selection procedure* overfits. `> 0.2` is a warning,
  `> 0.5` means picking the in-sample best is worse than random.
- **HAC t** (`qvalidate.hac_t`) — the naive t is inflated on overlapping/autocorrelated data.
- **MinTRL** (`qvalidate.min_track_record_length`) — periods needed before a Sharpe is credible.
- **FDR** (`qvalidate.benjamini_hochberg`) — when testing many hypotheses, nominal p-values are void.

### Build the null, don't assume it (`qmontecarlo.py`)
- **sign_flip_test** — the distribution-free test for the MEAN. A permutation of the order cannot
  change a sum, so use sign flips here, not reordering.
- **permutation_test** — for claims about timing/sequence (`sharpe`, `max_dd`), block-shuffled when
  the data is autocorrelated.
- **surrogate_paths** — simulate a no-edge process calibrated to the observed volatility
  (block bootstrap / iid / GBM / OU) and locate the real result in that distribution. Run at least
  two nulls of different character; if the verdict depends on which, that dependence is the result.
- **reality_check / spa_test** — the correct tests after screening many candidates against a
  benchmark. A per-strategy p-value is not.
- **stress_moments** — does the conclusion survive a 30% volatility understatement, doubled costs,
  the best period removed, a shrunk mean, an extra tail event? Publish the failing rows too.

### The best-of-N trap (deflated Sharpe / multiple testing)
If you searched `N` variants and picked the best, the winner's edge is inflated by
`≈ σ_total · √(2 ln N)`. With `σ_total = σ_period·√n`:

> **Search headroom** = `σ_period · √n · √(2 ln N)`

**What that quantity actually is, stated exactly** — because it is routinely used as if it
were a 5% test, and it is not. `√(2 ln N)` is the leading term of the asymptotic for the
max of N iid normals. It is **not** the mean of that max, and the significance level it
corresponds to **drifts with N**:

| N | 2 | 8 | 100 | 462 | 10 000 |
|---|---|---|---|---|---|
| `√(2 ln N)` | 1.18 | 2.04 | 3.04 | 3.50 | 4.29 |
| level it really screens at | 22.5% | **15.4%** | 11.3% | 10.1% | 8.5% |

So "beats the headroom" is a **weak screen**, not significance. For a stated level use
`qstats.search_threshold(N, alpha)` — exact, solving `1 − Φ(t)^N = α`:
`search_threshold(8) = 2.49` vs `√(2 ln 8) = 2.04`; `search_threshold(462) = 3.69` vs `3.50`.

**And the correction cuts both ways, so it is measured and not assumed.** Those thresholds
assume the N tries are *independent*. Real grids are correlated, so effective N is smaller
and the true threshold is **lower** — with mean ρ 0.3 across 462 alphas the effective
breadth collapses to ~3.3 and the 5% threshold falls to ~2.16. Use
`qportfolio.effective_breadth`; note it also **refuses** to call the number breadth when
ρ ≤ 0, and reports the cancellation signature instead.

Any result **below** the headroom is indistinguishable from luck. Report it:
> "Best of 60 variants, headroom ≈ ±2400pt; the winner beat baseline by 785pt → inside the noise."

Conversely — and this is the useful inverse — **if 0 of N variants beat the baseline while the
headroom was large, that is strong evidence the deficit is structural**, not sampling.

### Sample-size reality check (LLN, §8)
Needed periods for a given effect: `n ≈ (t_target · σ / effect)²`.
Run it *before* collecting, and say how long the wait is at the observed rate of new periods.
"15 sessions, ~1/weekday → 40-50 sessions is 5-7 weeks" is a decision-grade statement.

### The contraction trap (algebra masquerading as edge)
Any rule that **caps upside and downside** (trailing stop, partial TP, tighter stop, lower size)
is a **contraction**: `P&L_new ≈ α + β·P&L_base`, `0 < β < 1`.
A contraction **necessarily** beats the original on any window where the original is negative.
So "it wins on the recent bad stretch" carries **zero information** about merit. Test contractions
on the **full** sample, and never select one using a window chosen because the baseline lost there.

---

## 5b. CONDITIONAL EDGES — "it works better when X", and sizing up on it

**Fires whenever anyone proposes filtering or sizing a strategy by context**: confluences,
regime filters, "only trade when volatility is high", "size 3x on A+ setups", signal stacking.
This is the most common way a losing strategy is made to look profitable, and the failure is
invisible from the inside — the procedure returns a confident positive when nothing is there.

### The number that reframes everything

Searching for the best filter **manufactures** significance. Measured on a real case
(2026-08-03, mean-reversion on 16y NQ, 22 context variables, all singles and pairs = 968 rules):

| | t |
|---|---|
| best rule, real P&L | **+3.38** |
| best rule, **shuffled** P&L (median of 200) | **+3.02** |
| 95th percentile under the null | **+3.97** |

The discovery sat inside its own noise band. Frozen and applied to gold: **t = +0.15**, and the
second rule flipped sign to **−0.86**. Nothing was there.

> **A searched filter's t is uninterpretable until compared to the t the SAME SEARCH produces
> on data with no edge.** Not against 2.0. Against its own null.

Rough guide before you even run it: searching `N` rules puts the expected best-of-N near
`√(2 ln N)` σ — 968 rules ⇒ ~3.7σ, 7 200 rules ⇒ ~4.2σ. **Adding more candidate variables
raises the bar faster than it raises your chance of clearing it.** Prefer few pre-registered
hypotheses with an economic mechanism over a wide automatic sweep.

### The five gates — `scripts/qconfluence.py`, `confluence_gate()`

1. **Dose-response, not a binary filter.** Score each period by how many conditions align, and
   require the mean P&L to rise **monotonically** across buckets. A top-tercile-vs-rest split is
   the weakest possible evidence and the easiest for noise to fake; an ordered gradient is much
   harder. It is also the only shape that answers the sizing question — you need the **slope**,
   not one bucket's mean.
2. **Null calibration of the search** — shuffle the P&L, re-run the *entire* search, compare
   against that distribution's 95th percentile. `confluence_gate` does it literally: pass
   `build_score(masks, pnl) -> score`, the function that built your score, and gate 2 shuffles
   the P&L in blocks of 25 periods (`block`), rebuilds the score from each shuffle and
   recomputes the same slope t, two-sided. A score fixed before looking at this P&L takes
   `preregistered=True` instead. With neither, gate 2 is not run: the masks alone cannot say
   how a score was built from them. Measured on 200 datasets per shape, noise passes it
   between 3.5% and 7% of the time whatever the base mean, the rules' coverage, a two-valued
   P&L or persistence in the data. **A trend that score and P&L share still passes it** (188
   of 200): no shuffle can tell a shared trend from an edge, so gates 3 and 4 have to.
   `search_null` remains as a standalone tool; it tests each rule's mean P&L against zero, so
   centre the P&L before using it on its own.
3. **Out-of-sample** — select on the first part, apply frozen to the rest. No re-tuning.
4. **Independent replication** — same frozen rule, different instrument. In the measured case
   this is the gate that killed it, after gates 1-3 all passed.
5. **Concentration** — jackknife. An effect living in three periods is a story about three periods.

A gate **not executed** is not a gate passed: `confluence_gate` returns NON DECIDIBILE, never a
positive, when replication is missing, or when gate 2 has neither `build_score` nor
`preregistered=True`.

### If all five pass, then sizing — and only then

- Weights come from the **fitted slope**, capped (`sizing_from_slope`). Never from ratios of
  bucket means: those are unbounded and blow up exactly where the estimate is least reliable.
- Price the downside first with `sizing_cost(w) = mean(w)/√(mean(w²))` — the Sharpe multiplier
  **if the confluence is pure noise**. Measured on the real case: 1x/3x costs **13%** of Sharpe,
  0.5x/2x costs **19%**, 1x/5x costs **23%**. That is the premium paid for being wrong, and it
  is why the gates come first.

### Validate the validator

`qconfluence` ships with a **positive control**: inject a known effect and confirm the gates
pass it. This caught a real defect — `dose_response` used quantile edges, which on a discrete
confluence count (0-3) produce duplicate edges, empty buckets and `nan` means, so gate 1 failed
**even on a true injected effect**. A gate that rejects the true positive is not a strict gate,
it is a broken one. Always run both directions.

*(Recurring implementation trap, hit three times in one session: in NumPy `bool + bool` is
logical OR, not a count. `sum(mask_a + mask_b + mask_c)` silently yields 0/1 instead of 0-3,
collapsing the dose-response to two buckets. Cast with `.astype(int)` / `int(...)` first.)*

---

## 5c. "NO EDGE WORKS FOREVER" — the objection, made testable

**Fires whenever a full-sample negative is explained away by regime change**, or whenever
rotation between strategies is proposed. The objection is serious and usually right in spirit;
the inference drawn from it usually is not.

**Where it is right:** `t = Sharpe·√years` assumes stationarity, and edge decay is real
(McLean & Pontiff: ~58% of anomaly returns vanish post-publication). Averaging over 16 years
dilutes an edge alive for 4 of them.

**Where the argument breaks:** dilution lowers `t`; **it does not flip the sign**. 4 good years
+ 12 flat = a positive mean with a low t. A *negative* full-sample result needs years where the
strategy **loses**, not years where it stops working. So the objection makes three falsifiable
predictions — `scripts/qregime.py`, `regime_report()`:

| | prediction | function |
|---|---|---|
| **P1** | strong contiguous windows exist | `window_scan` |
| **P2** | **more** of them than chance gives | `window_scan` vs block null |
| **P3** | the regime is identifiable **in advance** | `regime_switch` |

**P2 is not optional.** Picking "the best window" is a search: 16 years scanned with 2-year
windows gives ~14 overlapping positions, and the max of 14 correlated windows beats the mean by
1.5-2σ with no edge at all. The null must be **block**-shuffled — an iid shuffle preserves
nothing and makes every series look regimed.

**P3 needs two baselines, not one.** A switch that is merely flat part of the time changes the
mean without skill. `regime_switch` always reports against always-on **and** a random switch at
the same duty cycle. Beating always-on but not random means the merit is in being out of the
market, not in the timing. (Observed in the built-in null control: `beats_always=True`,
`beats_random=False`.)

Measured, 2026-08-03, mean-reversion on 16y NQ: best 2-year window `t = +0.06`, `p = 0.983` —
98% of *shuffled* series produce a better best-window, and the trailing switch made it worse
(−128 → −333), worse than random too. Verdict: **never an edge**, not a diluted one.

`validate_power()` builds a synthetic edge confined to a quarter of the sample and confirms the
test fires on it (`p = 0.000`, switch beats both baselines). **Run it before trusting a
negative** — a regime test that cannot detect a regime-limited edge says nothing.

### And the practical consequence

Rotation does not lower the evidence bar — it raises it, because you must now also prove the
detector works ex ante. Where you cannot buy **years**, buy **breadth and frequency**:

```
t  = Sharpe × √years      ← time cannot be bought
IR = IC × √breadth        ← breadth can
```

A tiny edge over 200 instruments and 50 trades a day is decidable in months; the same edge on
one pair and one trade a day is never decidable. When a strategy is undecidable, the productive
move is usually to change the **source of observations**, not to weaken the standard.

---

## 6. CONTRARIAN PASS — mandate: destroy the finding

Run `references/contrarian-protocol.md` in full. Spawn an independent adversarial agent for
anything that moves capital, changes strategy, or goes to a third party. Its brief is to **refute**,
not to review. Report whether the finding survived, and what it cost.

Non-negotiable attacks (details in the protocol):
1. **Measurability** — is the rule below the data's resolution? (§2)
2. **Algebra** — is the "win" mechanical (contraction, leverage identity, denominator change)?
3. **Selection** — was the winner chosen post-hoc, on a window used for another criterion?
4. **Concentration** — does the result live in 1-3 periods? Jackknife it.
5. **Regime** — one regime = one observation. Which regime is missing (trend, crash, chop, high-rate)?
6. **Executability** — poll rate, trigger basis (mark vs last), min size, margin, borrow.
7. **Costs** — measured or assumed? adverse selection? funding on longer holds?
7b. **Fillability** — `qmicro.unfillable_share`. For any entry AT A LEVEL (breakout, stop,
   limit, a touch of a prior high/low): was the price already past it when the window opened?
   A London-range breakout showed t=+4.72 / +13.11 across two instruments, 16/17 positive
   years, jackknife 20/20 — and 36-44% of its "trades" were unfillable gap opens carrying
   **106% and 96%** of the profit. Real intraday breakouts: t −0.35 and +0.69. **No amount of
   statistical rigour catches a wrong fill assumption.**
8. **The missing number** — is the headline metric of the *recommended* variant actually reported,
   or was a different variant's number quoted? (This one catches real mistakes: always print the
   full metric set for the exact configuration being recommended.)

---

## 7. REQUIRED CHARTS

Tables hide shape; charts expose it. Produce these (see `charts.md`, `scripts/qcharts.py`),
and follow the `dataviz` skill for styling if a visual is being published:

1. **Equity curve** — baseline vs variant(s), same axes. Look for: is the gain a step (one period) or a slope?
2. **Underwater / drawdown curve** — depth AND duration. Two strategies with equal return and
   different time-under-water are not comparable.
3. **Per-period contribution bars** — sorted; instantly reveals concentration and which periods drive the delta.
4. **Parameter surface** — metric vs parameter. A **sharp isolated peak = overfit**; a broad plateau = robust.
   Never adopt a parameter at a lone spike; prefer the middle of a plateau.
5. **Return distribution / QQ vs normal** — where the tails actually are, and whether Sharpe/VaR
   (which assume normality) are lying (see `04-finance.md` §28, §30).

---

### Deliverable
Assemble the analysis with `qreport.Report` and hand over the HTML, not a terminal dump. The report
must contain the assumptions, the metric set for **every** variant (not just the recommended one),
the uncertainty, the charts, the verdict, and an explicit list of what the evidence does **not**
support. If a section would be empty, say so rather than omitting it silently.

## 8. VERDICT — say which of the three it is

1. **DECIDABLE, positive** — effect, band, t, why it survived the contrarian pass, what regime
   it has NOT seen. Give the sizing implication in currency and % of capital.
2. **DECIDABLE, negative** — the claim fails; state the mechanism that kills it (not just "not significant").
3. **NOT DECIDABLE** — say it plainly, then give the shortest path to decidability: what data,
   at what resolution, for how long, and the **pre-registered** criterion that will be applied to it.

Then the default action, explicitly: **change nothing** unless (1) is true. The burden of proof is
on the change, never on the status quo.

Close with the standard status line: **done · next · open risks**.

**Never present output as investment advice.** This protocol evaluates evidence; it does not
recommend positions to third parties, and backtested figures are hypothetical (see README disclaimer).

---

## 8b. COST IS A MEASUREMENT, NOT AN ASSUMPTION

An unverified cost model silently decides verdicts, and it is the assumption most likely to
be wrong by an order of magnitude. Before publishing any *net* figure:

1. Estimate the spread from the data (`qmicro.roll_spread`, `qmicro.corwin_schultz_spread`).
2. If short-horizon reversion is part of the thesis, run `qmicro.bounce_check`.
3. Express cost as a **fraction of the target**, and compute the frequency ceiling
   (`qmicro.max_trade_frequency`).
4. The first time realised P&L disagrees with the simulator, run
   `qmicro.decompose_slippage` — entry slip, exit slip and commission have three different
   fixes, and conflating them wastes months.
5. Any cost number not measured from real fills stays marked `⚠ assumed` in the output, and
   every conclusion that depends on it is stress-tested at 2×.

## 9. HOUSE RULES (adapt to your project)

- **Single source of numbers.** Every published figure must come from one designated script/query
  with its assumptions stated, or it doesn't get said. Check for drift by re-running it.
- **Bands, not points.** State the assumptions that move the band.
- **Base-rate discipline.** A previously retracted thesis is the base rate for the next one.
  Out-of-sample or nothing. Leverage that improves a backtest is a red flag until proven out-of-sample.
- **New strategies stay ring-fenced** from committed capital until a live, multi-regime track record exists.
- **Never quote an unvalidated edge** to investors, clients, or on a public page.
- **Never model taxes** unless that is explicitly the task.
- The analyst computes and verifies; a human signs and executes. No key handling, no order placement.
