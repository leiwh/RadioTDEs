Radio data of TDEs

A description is available at https://arxiv.org/abs/2603.06056

---

## k_CR calculator for radio TDEs

A single-file Python tool that selects the correct closure relation for the
circumnuclear-medium (CNM) density index `k_CR` and evaluates it from the
temporal index `α` and the electron spectral index `p`.

It implements every row of Tables 2-5 of the accompanying paper
([arXiv:2603.06056](https://arxiv.org/abs/2603.06056)): relativistic jets in
the coasting and Blandford-McKee (BM) phases, and non-relativistic outflows
in the coasting, Sedov-Taylor (ST) and deep-Newtonian (DN) phases.

Convention: `Fν ∝ t^(-α) ν^(-β)` and `n(R) ∝ R^(-k)`.

### Requirements

Python 3.8 or newer. Only the standard library is used.

### Usage

Interactive (Chinese by default, `-e` for English):

```bash
python3 kcr_calculator.py
python3 kcr_calculator.py -e
```

The script asks for `α`, `p`, the outflow type, the cooling regime, the
frequency ordering and the dynamical phase, then prints the selected table
row, the formula and `k_CR`.

Non-interactive:

```bash
python3 kcr_calculator.py -a -1.16 -p 2.28 --jetted --cooling slow --order amc --phase bm
python3 kcr_calculator.py -a 2.55 -p 3.04 --nonrel --order amc --phase coasting
```

An uncertainty may be appended as `+-` (or `-+`), e.g. `-a -1.16+-0.02`.
Use `--json` for machine-readable output.

### Options

| Option | Meaning |
| --- | --- |
| `-a, --alpha` | temporal index `α` (uncertainty optional) |
| `-p, --p` | electron spectral index `p` (uncertainty optional) |
| `--jetted` / `--nonrel` | outflow type: `β = -2` or `β = (p-1)/2` |
| `--cooling {fast,slow}` | jetted only: `νa<νc<νm` or slow cooling |
| `--order {acm,amc,mac}` | frequency ordering |
| `--phase {coasting,bm,st,dn,both}` | dynamical phase |
| `-e, --english` | English output (default: Chinese) |
| `--json` | JSON output |

### Self-checks

```bash
python3 kcr_calculator.py --selftest      # inverts alpha(k) back to k for every row
python3 kcr_calculator.py --verify-paper  # reproduces the 26 k_CR values of Table 6
python3 kcr_calculator.py --list-tables   # prints all formulas of Tables 2-5
```

### Notes

- `β` is set by the outflow type: `β = -2` for jetted TDEs (low-frequency,
  self-absorbed plateau) and `β = (p-1)/2` for non-relativistic outflows.
- For jetted TDEs in slow cooling, the coasting rows carry no `k_CR`
  (`α = -2` independent of `k`); the script reports this instead of a number.
- Table 5 (deep-Newtonian) has two branches, `2 < p < 3` and `p >~ 3`.

### Reference

If you use this script, please cite [arXiv:2603.06056](https://arxiv.org/abs/2603.06056).
