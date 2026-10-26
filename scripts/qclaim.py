#!/usr/bin/env python3
"""qclaim — le affermazioni si GENERANO dai numeri, non si scrivono a mano.

PERCHE' ESISTE. Nel luglio-agosto 2026 lo stesso errore e' arrivato all'output cinque volte,
sempre nella stessa forma: **una frase che contraddice una tabella stampata due righe sopra**.

    "11.71 resta sotto 8.37"          → 11.71 > 8.37, era stampato sopra
    "salva il 60-70%"                 → misurato 25%, era stampato sopra
    "3.50 e' la media del rumore"     → la media era 3.013, era stampato sopra
    "aggiornata la casella del gate"  → grep -c diceva 0, era stampato sopra
    "il costo/range decide"           → metrica sbagliata per uno strumento binario

In OGNI caso il numero calcolato era giusto e la frase sbagliata, e la smentita era gia'
visibile nell'output. Non e' un problema di conoscenza: e' che ho RICORDATO invece di
GUARDARE, a due righe di distanza.

Regole in prosa scritte per impedirlo: R10, R11, G8, G9. **Violate tutte entro 48 ore.**
Guardie eseguibili scritte nello stesso periodo (`reachable_threshold`, `plausibility_flag`,
`borrowed_parameter`, il controllo `tickv == tick*mult`): **hanno preso ogni bug.**

La differenza non e' la qualita' della regola. E' che una funzione scatta da sola e una
frase no. Quindi qui non c'e' una regola nuova: c'e' il modo di non poter sbagliare.

USO MINIMO, prima di consegnare qualunque affermazione comparativa:

    from qclaim import cmp_, rel, audit
    print(cmp_(breadth_new, breadth_old, "breadth 34 strumenti", "breadth 6 gambe"))
    print(rel(t_stat, ">=", 2.0, "t del risultato principale"))
    audit(draft_text, {"media": 3.013, "soglia5": 3.693})   # numeri non tracciabili → FAIL
"""
import math
import re

__all__ = ["cmp_", "rel", "pct_change", "audit", "verdict", "ClaimError"]


class ClaimError(AssertionError):
    """Un'affermazione non regge ai numeri che dovrebbe descrivere."""


# ---------------------------------------------------------------- 65. confronto generato


def cmp_(a, b, la="A", lb="B", unit="", fmt="{:.4g}", strict=True):
    """La frase di confronto, GENERATA dai valori. Impossibile invertirla per distrazione.

    Il caso reale: ho scritto "11.71 resta sotto 8.37" due righe sotto una tabella che
    mostrava 11.71 e 8.37. Con questa funzione la frase esce da sola e dice "sopra".
    """
    a = float(a); b = float(b)
    fa, fb = fmt.format(a) + unit, fmt.format(b) + unit
    if math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15):
        word = "UGUALE a"
    elif a > b:
        word = "SOPRA"
    else:
        word = "SOTTO"
    delta = a - b
    ratio = (a / b) if b not in (0,) else float("nan")
    txt = f"{la} = {fa} è {word} {lb} = {fb}  (Δ {fmt.format(delta)}{unit}"
    if b != 0 and (a * b) > 0:
        txt += f" · {ratio:.3g}×"
    txt += ")"
    return txt


def rel(value, op, threshold, label="valore", unit="", fmt="{:.4g}", raise_on_fail=False):
    """Verifica una relazione e restituisce la frase. `op` in >, >=, <, <=, ==, !=.

    Serve per le affermazioni di soglia — "t supera 2", "il costo sta sotto il target" —
    che sono la seconda fonte di contraddizioni fra prosa e tabella.
    """
    v = float(value); t = float(threshold)
    ops = {">": v > t, ">=": v >= t, "<": v < t, "<=": v <= t,
           "==": math.isclose(v, t, rel_tol=1e-12, abs_tol=1e-15),
           "!=": not math.isclose(v, t, rel_tol=1e-12, abs_tol=1e-15)}
    if op not in ops:
        raise ValueError(f"operatore non valido: {op}")
    ok = ops[op]
    txt = (f"{label} = {fmt.format(v)}{unit} {op} {fmt.format(t)}{unit} → "
           f"{'VERO' if ok else 'FALSO'}")
    if not ok and raise_on_fail:
        raise ClaimError(txt)
    return txt if ok else txt + "  ⚠ l'affermazione NON regge"


def pct_change(new, old, label="variazione", fmt="{:.1f}"):
    """Variazione percentuale col segno giusto. Evita gli 'aumenta del 25%' su un calo."""
    new = float(new); old = float(old)
    if old == 0:
        return f"{label}: base zero, variazione percentuale non definita"
    ch = 100.0 * (new - old) / abs(old)
    word = "aumenta" if ch > 0 else ("cala" if ch < 0 else "invariato")
    return f"{label}: {word} del {fmt.format(abs(ch))}% ({fmt.format(old)} → {fmt.format(new)})"


# ---------------------------------------------------------------- 66. tracciabilità


_NUM = re.compile(r"[-+]?\d{1,3}(?:[ .,]\d{3})*(?:[.,]\d+)?|[-+]?\d*[.,]\d+|[-+]?\d+")


def _parse(tok):
    t = tok.strip().replace(" ", "")
    if t.count(",") and t.count("."):
        t = t.replace(",", "") if t.rfind(".") > t.rfind(",") else t.replace(".", "").replace(",", ".")
    elif t.count(","):
        t = t.replace(",", ".") if len(t.split(",")[-1]) != 3 else t.replace(",", "")
    try:
        return float(t)
    except ValueError:
        return None


def audit(text, values, rel_tol=0.02, ignore=(), min_abs=0.0):
    """Ogni numero nella prosa è TRACCIABILE a una quantità calcolata?

    Il caso reale: ho scritto "salva il 60-70%" dove la misura diceva 25%. Il 60 e il 70
    non venivano da nessun calcolo — li avevo ricordati. Questa funzione li avrebbe marcati.

    `values`: dict {etichetta: valore} di tutto ciò che hai davvero calcolato.
    Ritorna dict con `ok`, `untraceable` (numeri senza origine) e `matched`.
    Numeri banali (0, 1, 2, 100, anni, indici) si escludono con `ignore` o `min_abs`.
    """
    vals = [float(v) for v in values.values() if isinstance(v, (int, float))]
    ign = set(float(x) for x in ignore)
    found, untraced = [], []
    for tok in _NUM.findall(text):
        x = _parse(tok)
        if x is None or x in ign or abs(x) < min_abs:
            continue
        hit = None
        for lab, v in values.items():
            if not isinstance(v, (int, float)):
                continue
            v = float(v)
            if math.isclose(x, v, rel_tol=rel_tol, abs_tol=1e-12) or \
               math.isclose(x, 100 * v, rel_tol=rel_tol, abs_tol=1e-12) or \
               math.isclose(x, v / 100, rel_tol=rel_tol, abs_tol=1e-12):
                hit = lab; break
        (found if hit else untraced).append((x, hit))
    return {"ok": not untraced, "n_checked": len(found) + len(untraced),
            "matched": found, "untraceable": [x for x, _ in untraced],
            "note": ("ogni numero risale a una quantità calcolata" if not untraced else
                     f"⚠ {len(untraced)} numeri NON tracciabili: {[x for x,_ in untraced]} — "
                     "ricordati o inventati? Ricalcolali o toglili.")}


# ---------------------------------------------------------------- 67. verdetto coerente


def verdict(metrics, rules, default="NON DECIDIBILE"):
    """Il verdetto si DERIVA dai criteri dichiarati, non si sceglie a fine analisi.

    `rules`: lista di (nome, callable(metrics)->bool, esito_se_vero), valutate in ordine.
    Serve contro l'ultima forma dell'errore: dichiarare un criterio a inizio lavoro
    ("falsificata se t < 1.6") e poi consegnare un verdetto che non lo applica.
    """
    trace = []
    for name, test, out in rules:
        try:
            r = bool(test(metrics))
        except Exception as e:
            r = False
            name = f"{name} (errore: {type(e).__name__})"
        trace.append((name, r))
        if r:
            return {"verdict": out, "matched": name, "trace": trace}
    return {"verdict": default, "matched": None, "trace": trace}
