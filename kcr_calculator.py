#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kcr_calculator.py — Radio TDE 的 CNM 密度斜率 k_CR 计算器（closure relations）
==============================================================================

依据
----
RadioTDECNM.pdf / ms.tex, Sec. 3.1 "Closure relations (CRs)" 的 Table 2-5：

    Table 2 : 相对论喷流 coasting 相（非相对论 coasting 相与之共享同一套 CR）
    Table 3 : BM（Blandford–McKee）减速相（相对论喷流）
    Table 4 : ST（Sedov–Taylor）减速相（非相对论外流）
    Table 5 : Deep-Newtonian（DN）减速相（非相对论外流）

约定（论文 Eq. 1）
------------------
    F_ν ∝ t^{-α} ν^{-β}

    α : 时间指数（流量上升时 α<0，衰减时 α>0）
    β : 谱指数
    k : CNM 密度剖面指数，n(R) ∝ R^{-k}
    p : 电子能谱指数，N(γ) ∝ γ^{-p}

    注意：论文中 α 带负号（F ∝ t^{-α}），所以 "α<1" 判别与常见文献里
    "F ∝ t^{α}" 的定义差一个负号，本代码严格采用论文定义。

选择规则（与需求一致）
----------------------
0) 先输入 α 和 p。
1) jetted TDE（相对论喷流）→ β = -2：
   1a) fast cooling（νa<νc<νm）：问 coasting 还是 BM，
       → 取该 block 中 β=-2 的行（ν<νa）；
   1b) slow cooling：先问 νa<νm<νc 还是 νm<νa<νc，再问 coasting 还是 BM，
       → 取对应 block 中 β=-2 的行（ν<νa 或 ν<νm）。
2) 非相对论外流 → β = (p-1)/2：
   先问 νa<νm<νc 还是 νm<νa<νc；
   若 α<1 → 用 coasting 动力学的 k_CR 表达式；
   若 α>=1 → 选择减速相 ST 或 Deep-Newtonian（也可两者都算）。

用法
----
交互式（按上面的问题顺序逐项询问）::

    python3 kcr_calculator.py

非交互式（脚本/批处理）::

    python3 kcr_calculator.py -a -1.16 -p 2.28 --jetted --cooling slow --order amc --phase bm
    python3 kcr_calculator.py -a 2.55  -p 3.04 --nonrel --order amc --phase coasting
    python3 kcr_calculator.py -a -1.16+-0.02 -p 2.28 --jetted --cooling slow --order amc --phase bm
    python3 kcr_calculator.py -a -1.16 -p 2.28 --jetted --cooling slow --order amc --phase bm --json

误差写法：``值+-误差`` 或 ``值-+误差``，例如 ``-1.16+-0.02``、``-1.16-+0.02``；
误差只用于数值传递（δk），不影响 k_CR 的中心值。

自检 / 工具::

    python3 kcr_calculator.py --selftest      # 用 α(k) 反解校验 Table 2-5 每一行公式
    python3 kcr_calculator.py --verify-paper  # 用论文 Table 6 的 26 个源复核 k_CR
    python3 kcr_calculator.py --list-tables   # 打印 Table 2-5 的全部公式

输出语言
--------
默认中文；加 ``-e`` / ``--english`` 显示英文（也适用于 ``--help`` 与交互式）::

    python3 kcr_calculator.py -e
    python3 kcr_calculator.py -e -a -1.16 -p 2.28 --jetted --cooling slow --order amc --phase bm

(Output is Chinese by default; add ``-e`` / ``--english`` for English.)

依赖：仅标准库（Python 3.8+）。
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

# --------------------------------------------------------------------------- #
# 0. 国际化：默认中文，-e / --english 切换英文
# --------------------------------------------------------------------------- #
_LANG = "zh"

STRINGS: Dict[str, Dict[str, str]] = {
    "zh": {
        # ---- report ----
        "rep_banner": " Radio TDE closure-relation calculator  —  k_CR from α and β",
        "rep_input": "输入 :  α = {a}      p = {b}",
        "rep_path": "路径 :  {p}",
        "rep_selected": "选中 :  {w}",
        "rep_beta": "        β = {b}   （该行 β 标签: {lab}）",
        "rep_alphak": "α(k)      = {f}",
        "rep_kcrformula": "k_CR(α,β) = {f}",
        "rep_k": "k_CR = {k}{err}",
        "rep_selfcheck": "自洽  : α(k_CR) = {m}  (输入 {a}, 残差 {r})",
        "rep_undet": "k_CR = —   （无法确定）",
        "rep_checks": "检查 :",
        "rep_candidates": "同一 β = {b} 下其它动力学阶段/排序的 k_CR（供对比）:",
        "rep_mark": "  ← 采用",
        "rep_foot_jetted": "  （Table 4/5 的 ST、Deep-Newtonian 仅适用于非相对论外流，故未列出）",
        "rep_foot_nonrel": "  （Table 3 的 BM 仅适用于相对论喷流，故未列出）",
        "tip_close": "        提示: 该行预测 α = {p}，输入 α = {a} 与之接近，说明该参数组合自洽；但 k 仍不可定。",
        "tip_off": "        提示: 该行预测 α = {p}，与输入 α = {a} 不符，该分支很可能不适用（可试其它阶段/冷却类型）。",
        "tip_singular": "        提示: 该行 α(k) 为 k=0→{a0}，k=3→{a3}；α={a} 使公式分母为零，请检查 α 取值。",
        # ---- warnings ----
        "w_p_le2": "p={p} ≤ 2：电子能谱过软，(p-1)/2 及 Deep-Newtonian 分支假设失效",
        "w_p_gt4": "p={p} > 4：超出论文 MCMC 先验范围 p∈[2,4]",
        "w_p_near3": "p={p} 接近 3：Table 5 的 2<p<3 与 p≳3 两分支在此交界，结果对 p 敏感",
        "w_k_neg": "k_CR={k} < 0：密度随半径增大（反转剖面），物理上通常不合理",
        "w_k_gt3": "k_CR={k} > 3：超出常见范围 k∈[0,3]（0=均匀，2=Bondi，3=星风）",
        "w_k_ok": "k_CR={k} ∈ [0,3]，处于常见范围",
        "w_alpha_extrap": "α={a} 超出 k∈[0,3] 在该行对应的 α 区间 [{lo}, {hi}] → k 为外插，谨慎使用",
        "w_kcr_none": "Table 中该行 k_CR 为 '−'：α = {lab} 与 k 无关，无法由闭合关系定出 k",
        "w_jetted_coasting": "论文 Sec 3.2.2 指出相对论 TDE 多处于 BM 减速相；coasting 结果仅供对比",
        "w_nonrel_alt": "注意: coasting 分支也给出 k={k}∈[0,3]；α≥1 本身不能排除 coasting，需结合动力学先验判断",
        "w_dn_amc": "Table 5 按 νm<νa<νc 推导；νa<νm<νc 时其 νm<ν<νc 段的 α(k) 与 νa<ν<νc 段相同，故该公式通用",
        # ---- messages / errors ----
        "msg_kcr_none": "Table 中该行 k_CR 列为 '−'：α = {lab} 是常数、与 k 无关，故无法由闭合关系定出 k",
        "msg_zero_den": "分母为零：该 α 使公式奇异",
        "msg_not_finite": "结果为非有限值（公式在该 α 处奇异）",
        "err_outflow": "未知外流类型: {v}",
        "err_phase": "未知动力学阶段: {v}",
        "err_ordering": "未知频率排序: {v}（应为 ACM/AMC/MAC）",
        "err_row_not_found": "在 {tbl} / {ord} 中找不到 β={b} 的行（p={p}）。请检查外流类型/冷却类型/α 范围。",
        "err_row_ambiguous": "{tbl} 中 β={b} 匹配到 {n} 行，无法唯一确定：{segs}",
        "where_row": "行 \"{seg}\"",
        # ---- labels ----
        "ph_coasting": "coasting（Table 2）",
        "ph_BM": "BM 减速（Table 3）",
        "ph_ST": "ST 减速（Table 4）",
        "ph_DN": "Deep-Newtonian 减速（Table 5）",
        "path_jetted": "jetted TDE",
        "path_nonrel": "非相对论外流",
        "path_fast": "fast cooling",
        "path_slow": "slow cooling",
        # ---- interactive ----
        # ---- interactive（精简版）----
        "it_title": "Radio TDE k_CR calculator  —  RadioTDECNM Tables 2-5",
        "it_err_fmt": "α, p: 可带误差，如 -1.16+-0.02",
        "it_alpha": "α = ",
        "it_p": "p = ",
        "it_ask_outflow": "1) 外流",
        "it_opt_jetted": "jetted",
        "it_opt_nonrel": "非相对论",
        "it_ask_cooling": "2) 冷却",
        "it_opt_fast": "fast",
        "it_opt_slow": "slow",
        "it_ask_order": "{n}) 排序",
        "it_opt_amc": "νa<νm<νc",
        "it_opt_mac": "νm<νa<νc",
        "it_ask_phase": "{n}) 阶段",
        "it_opt_coasting": "coasting (T2)",
        "it_opt_bm": "BM (T3)",
        "it_alpha_lt1": "α = {a} < 1 → coasting (T2)",
        "it_ask_decel": "3) 减速相",
        "it_opt_st": "ST (T4)",
        "it_opt_dn": "DN (T5)",
        "it_opt_both": "两个都算",
        "it_bad_choice": "  无效，请重选",
        "it_bad_number": "  无法解析，如 -1.16+-0.02",
        # ---- selftest / verify ----
        "st_summary": "公式自检: {p} 通过, {f} 失败, {s} 跳过（奇异点）",
        "st_ok": "  ✓ 全部通过：Table 2-5 的 α(k) 与 k_CR(α,β) 互相一致",
        "st_bad": "  ✗ 存在不一致，请检查 Registry",
        "st_fail_none": "  FAIL {w}: k_CR='−' 但 α 依赖 k",
        "st_fail_inv": "  FAIL {w}  p={p} k={k}: α={a} → k_back={kb}",
        "vp_header": "用论文 Table 6（26 个源）复核 k_CR 公式与 label 映射",
        "vp_legend": "label: a1/a2 = BM 减速(a1:νa<νm<νc, a2:νm<νa<νc); b1/b2 = coasting(b1:νa<νm<νc, b2:νm<νa<νc)",
        "vp_cols": "k_表={kt}  k_算={kc}  Δ={d}",
        "vp_summary": "结果: {n}/{N} 通过（容差 0.02，属于论文表格舍入）",
        # ---- CLI ----
        "ap_desc": "Radio TDE closure-relation calculator: 由 α, p 选取并计算 k_CR（Table 2-5）",
        "ap_epilog": ("示例:\n"
                      "  python3 kcr_calculator.py                       # 交互式（中文）\n"
                      "  python3 kcr_calculator.py -e                    # 交互式（英文）\n"
                      "  python3 kcr_calculator.py -a -1.16 -p 2.28 --jetted --cooling slow --order amc --phase bm\n"
                      "  python3 kcr_calculator.py -a 2.55  -p 3.04 --nonrel --order amc --phase coasting\n"
                      "  python3 kcr_calculator.py --selftest | --verify-paper | --list-tables"),
        "ap_alpha": "时间指数 α（F∝t^-α），可写 -1.16 或 -1.16+-0.02",
        "ap_p": "电子能谱指数 p，可写 2.5 或 2.5+-0.1",
        "ap_jetted": "jetted TDE（β=-2）",
        "ap_nonrel": "非相对论外流（β=(p-1)/2）",
        "ap_cooling": "jetted 时：fast cooling / slow cooling",
        "ap_order": "频率排序（ACM=νa<νc<νm; AMC=νa<νm<νc; MAC=νm<νa<νc）",
        "ap_phase": "动力学阶段",
        "ap_json": "以 JSON 输出",
        "ap_selftest": "公式自检（α(k) 反解）",
        "ap_verify": "用论文 Table 6 的 26 个源复核",
        "ap_list": "打印 Table 2-5 全部公式",
        "ap_english": "显示英文（默认中文）",
        "cli_need_outflow": "非交互模式必须给出 --jetted 或 --nonrel",
        "cli_need_cooling": "jetted 模式必须给出 --cooling fast|slow",
        "cli_need_phase": "jetted 模式必须给出 --phase coasting|bm",
        "cli_note_fast": "提示: fast cooling 对应 νa<νc<νm，忽略 --order",
        "cli_need_order_slow": "jetted + slow cooling 需给出 --order amc|mac",
        "cli_need_order_nonrel": "非相对论模式需给出 --order amc|mac",
        "cli_need_phase_alpha1": "α≥1 时需给出 --phase st|dn|both",
        "cli_need_alpha_p": "非交互模式需同时给出 -a/--alpha 与 -p/--p",
    },
    "en": {
        # ---- report ----
        "rep_banner": " Radio TDE closure-relation calculator  —  k_CR from α and β",
        "rep_input": "Input   :  α = {a}      p = {b}",
        "rep_path": "Path    :  {p}",
        "rep_selected": "Selected:  {w}",
        "rep_beta": "        β = {b}   (row β label: {lab})",
        "rep_alphak": "α(k)      = {f}",
        "rep_kcrformula": "k_CR(α,β) = {f}",
        "rep_k": "k_CR = {k}{err}",
        "rep_selfcheck": "Check : α(k_CR) = {m}  (input {a}, residual {r})",
        "rep_undet": "k_CR = —   (undetermined)",
        "rep_checks": "Checks:",
        "rep_candidates": "k_CR for other phases/orderings at the same β = {b} (for comparison):",
        "rep_mark": "  <- selected",
        "rep_foot_jetted": "  (ST and Deep-Newtonian in Tables 4-5 apply to non-relativistic outflows only; not listed)",
        "rep_foot_nonrel": "  (BM in Table 3 applies to relativistic jets only; not listed)",
        "tip_close": "        Note: this row predicts α = {p}, close to the input α = {a}; the parameter set is self-consistent, but k remains undetermined.",
        "tip_off": "        Note: this row predicts α = {p}, which does not match the input α = {a}; this branch is likely inapplicable (try another phase / cooling type).",
        "tip_singular": "        Note: this row has α(k) = {a0} at k=0 and {a3} at k=3; α = {a} makes the denominator vanish — please check α.",
        # ---- warnings ----
        "w_p_le2": "p={p} <= 2: electron spectrum too soft; the (p-1)/2 and Deep-Newtonian branches break down",
        "w_p_gt4": "p={p} > 4: outside the paper's MCMC prior p in [2, 4]",
        "w_p_near3": "p={p} near 3: the 2<p<3 and p>~3 branches of Table 5 meet here; the result is sensitive to p",
        "w_k_neg": "k_CR={k} < 0: density increases with radius (inverted profile); usually unphysical",
        "w_k_gt3": "k_CR={k} > 3: outside the usual range k in [0, 3] (0=uniform, 2=Bondi, 3=wind)",
        "w_k_ok": "k_CR={k} in [0, 3]: within the usual range",
        "w_alpha_extrap": "α={a} lies outside the α range [{lo}, {hi}] spanned by k in [0, 3] for this row -> k is an extrapolation; use with caution",
        "w_kcr_none": "This row has no k_CR in the table ('-'): α = {lab} is independent of k, so k cannot be determined from the closure relation",
        "w_jetted_coasting": "Sec. 3.2.2 of the paper notes that relativistic TDEs are mostly in the BM deceleration phase; the coasting result is for comparison only",
        "w_nonrel_alt": "Note: the coasting branch also gives k={k} in [0, 3]; α >= 1 alone does not exclude coasting, so a dynamical prior is needed",
        "w_dn_amc": "Table 5 is derived for νm<νa<νc; for νa<νm<νc the νm<ν<νc segment has the same α(k) as the νa<ν<νc segment, so the formula applies",
        # ---- messages / errors ----
        "msg_kcr_none": "This row's k_CR column is '-': α = {lab} is constant and independent of k, so k cannot be determined from the closure relation",
        "msg_zero_den": "Denominator is zero: the formula is singular at this α",
        "msg_not_finite": "Non-finite result (the formula is singular at this α)",
        "err_outflow": "Unknown outflow type: {v}",
        "err_phase": "Unknown dynamical phase: {v}",
        "err_ordering": "Unknown spectral ordering: {v} (expected ACM/AMC/MAC)",
        "err_row_not_found": "No row with β={b} found in {tbl} / {ord} (p={p}). Check the outflow type / cooling type / α range.",
        "err_row_ambiguous": "{tbl}: β={b} matches {n} rows; cannot choose uniquely: {segs}",
        "where_row": "row \"{seg}\"",
        # ---- labels ----
        "ph_coasting": "coasting (Table 2)",
        "ph_BM": "BM deceleration (Table 3)",
        "ph_ST": "ST deceleration (Table 4)",
        "ph_DN": "Deep-Newtonian deceleration (Table 5)",
        "path_jetted": "jetted TDE",
        "path_nonrel": "non-relativistic outflow",
        "path_fast": "fast cooling",
        "path_slow": "slow cooling",
        # ---- interactive ----
        # ---- interactive (compact) ----
        "it_title": "Radio TDE k_CR calculator  —  RadioTDECNM Tables 2-5",
        "it_err_fmt": "α, p: error optional, e.g. -1.16+-0.02",
        "it_alpha": "α = ",
        "it_p": "p = ",
        "it_ask_outflow": "1) outflow",
        "it_opt_jetted": "jetted",
        "it_opt_nonrel": "non-rel",
        "it_ask_cooling": "2) cooling",
        "it_opt_fast": "fast",
        "it_opt_slow": "slow",
        "it_ask_order": "{n}) ordering",
        "it_opt_amc": "νa<νm<νc",
        "it_opt_mac": "νm<νa<νc",
        "it_ask_phase": "{n}) phase",
        "it_opt_coasting": "coasting (T2)",
        "it_opt_bm": "BM (T3)",
        "it_alpha_lt1": "α = {a} < 1 -> coasting (T2)",
        "it_ask_decel": "3) deceleration",
        "it_opt_st": "ST (T4)",
        "it_opt_dn": "DN (T5)",
        "it_opt_both": "both",
        "it_bad_choice": "  invalid, try again",
        "it_bad_number": "  cannot parse; e.g. -1.16+-0.02",
        # ---- selftest / verify ----
        "st_summary": "Formula self-test: {p} passed, {f} failed, {s} skipped (singular)",
        "st_ok": "  ✓ All passed: α(k) and k_CR(α,β) are mutually consistent for every row of Tables 2-5",
        "st_bad": "  ✗ Inconsistent rows found; check the registry",
        "st_fail_none": "  FAIL {w}: k_CR='-' but α depends on k",
        "st_fail_inv": "  FAIL {w}  p={p} k={k}: α={a} → k_back={kb}",
        "vp_header": "Cross-check k_CR formulas and label mapping against Table 6 (26 sources)",
        "vp_legend": "label: a1/a2 = BM deceleration (a1: νa<νm<νc, a2: νm<νa<νc); b1/b2 = coasting (b1: νa<νm<νc, b2: νm<νa<νc)",
        "vp_cols": "k_tab={kt}  k_cal={kc}  Δ={d}",
        "vp_summary": "Result: {n}/{N} passed (tolerance 0.02, paper-table rounding)",
        # ---- CLI ----
        "ap_desc": "Radio TDE closure-relation calculator: select and evaluate k_CR from α and p (Tables 2-5)",
        "ap_epilog": ("Examples:\n"
                      "  python3 kcr_calculator.py                       # interactive (Chinese)\n"
                      "  python3 kcr_calculator.py -e                    # interactive (English)\n"
                      "  python3 kcr_calculator.py -a -1.16 -p 2.28 --jetted --cooling slow --order amc --phase bm\n"
                      "  python3 kcr_calculator.py -a 2.55  -p 3.04 --nonrel --order amc --phase coasting\n"
                      "  python3 kcr_calculator.py --selftest | --verify-paper | --list-tables"),
        "ap_alpha": "temporal index α (F ∝ t^-α); write -1.16 or -1.16+-0.02",
        "ap_p": "electron spectral index p; write 2.5 or 2.5+-0.1",
        "ap_jetted": "jetted TDE (β=-2)",
        "ap_nonrel": "non-relativistic outflow (β=(p-1)/2)",
        "ap_cooling": "for jetted: fast cooling / slow cooling",
        "ap_order": "frequency ordering (ACM=νa<νc<νm; AMC=νa<νm<νc; MAC=νm<νa<νc)",
        "ap_phase": "dynamical phase",
        "ap_json": "output JSON",
        "ap_selftest": "formula self-test (invert α(k))",
        "ap_verify": "cross-check with the 26 sources of Table 6",
        "ap_list": "print all formulas of Tables 2-5",
        "ap_english": "show English (default: Chinese)",
        "cli_need_outflow": "non-interactive mode requires --jetted or --nonrel",
        "cli_need_cooling": "jetted mode requires --cooling fast|slow",
        "cli_need_phase": "jetted mode requires --phase coasting|bm",
        "cli_note_fast": "Note: fast cooling corresponds to νa<νc<νm; --order ignored",
        "cli_need_order_slow": "jetted + slow cooling requires --order amc|mac",
        "cli_need_order_nonrel": "non-relativistic mode requires --order amc|mac",
        "cli_need_phase_alpha1": "α >= 1 requires --phase st|dn|both",
        "cli_need_alpha_p": "non-interactive mode requires both -a/--alpha and -p/--p",
    },
}


def set_language(code: str) -> None:
    """设置输出语言：'zh'（默认）或 'en'。"""
    global _LANG
    _LANG = code if code in STRINGS else "zh"


def get_language() -> str:
    return _LANG


def tr(key: str, **kw) -> str:
    """取当前语言的文案；缺失时回退中文。"""
    s = STRINGS.get(_LANG, STRINGS["zh"]).get(key) or STRINGS["zh"].get(key, key)
    return s.format(**kw) if kw else s


# --------------------------------------------------------------------------- #
# 1. 标签与常数
# --------------------------------------------------------------------------- #
ACM, AMC, MAC = "ACM", "AMC", "MAC"          # 频率排序
COASTING, BM, ST, DN = "coasting", "BM", "ST", "DN"

ORDER_LABEL: Dict[str, str] = {
    ACM: "νa<νc<νm  (fast cooling)",
    AMC: "νa<νm<νc  (slow cooling)",
    MAC: "νm<νa<νc  (slow cooling)",
}
ORDER_LABEL_PLAIN: Dict[str, str] = {
    ACM: "νa<νc<νm",
    AMC: "νa<νm<νc",
    MAC: "νm<νa<νc",
}
_PHASE_KEY: Dict[str, str] = {
    COASTING: "ph_coasting", BM: "ph_BM", ST: "ph_ST", DN: "ph_DN",
}


def phase_label(phase: str) -> str:
    """动力学阶段的（本地化）显示名。"""
    return tr(_PHASE_KEY[phase])
TABLE_NAME: Dict[int, str] = {
    2: "Table 2 (coasting)",
    3: "Table 3 (BM deceleration)",
    4: "Table 4 (ST deceleration)",
    5: "Table 5 (deep-Newtonian)",
}
BETA_TOL = 1e-6


# --------------------------------------------------------------------------- #
# 1. Table 2-5 的数据结构
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Row:
    """Table 2-5 中的一行。

    beta_fn(p)      : 该行的 β（-2、-1/3、1/2、-5/2、p/2、(p-1)/2）
    alpha_of_k(k,p) : 该行的 α(k)（论文表格第 3 列）
    kcr(a,b,p)      : 该行的 k_CR(α,β)（论文表格第 5 列）；None 表示表中为 "−"
    """
    table: int
    phase: str
    ordering: str
    block: int                       # 表内第几个 block（1 起）
    segment: str                     # 频率区间，如 "ν<νa"
    beta_fn: Callable[[float], float]
    beta_label: str
    alpha_of_k: Callable[[float, float], float]
    alpha_label: str
    kcr: Optional[Callable[[float, float, float], float]]
    kcr_label: str
    p_range: Optional[Tuple[float, float]] = None
    note: str = ""

    # -- 便捷方法 ---------------------------------------------------------
    def beta(self, p: float) -> float:
        return self.beta_fn(p)

    def in_p_range(self, p: float) -> bool:
        if self.p_range is None:
            return True
        lo, hi = self.p_range
        return lo <= p < hi

    def match_beta(self, beta: float, p: float, tol: float = BETA_TOL) -> bool:
        return abs(self.beta(p) - beta) < tol

    def where(self) -> str:
        """人类可读的定位串，方便回查论文表格。"""
        return (f"{TABLE_NAME[self.table]}, block {self.block} "
                f"({ORDER_LABEL_PLAIN.get(self.ordering, self.ordering)}), "
                + tr("where_row", seg=self.segment))


def _c(v: float) -> Callable[[float], float]:
    """常数 β。"""
    return lambda p, _v=v: _v


def _p_half(p: float) -> float:
    return (p - 1.0) / 2.0


def _p_over_2(p: float) -> float:
    return p / 2.0


def R(table: int, phase: str, ordering: str, block: int, seg: str,
      bfn: Callable[[float], float], blab: str,
      afn: Callable[[float, float], float], alab: str,
      kfn: Optional[Callable[[float, float, float], float]], klab: str,
      prange: Optional[Tuple[float, float]] = None, note: str = "") -> Row:
    return Row(table, phase, ordering, block, seg, bfn, blab, afn, alab, kfn, klab, prange, note)


# ============================ Table 2: coasting ============================ #
_T2_AMC_NOTE = "α=-2 与 k 无关 → 该行无法定 k"
ROWS: List[Row] = [
    # ---- Table 2, block 1: νa<νc<νm (fast cooling) ----
    R(2, COASTING, ACM, 1, "ν<νa", _c(-2.0), "-2",
      lambda k, p: -(1.0 + k), "-(1+k)",
      lambda a, b, p: (2 * a - b) / b, "(2α-β)/β"),
    R(2, COASTING, ACM, 1, "νa<ν<νc", _c(-1.0 / 3.0), "-1/3",
      lambda k, p: -(11.0 - 6.0 * k) / 3.0, "-(11-6k)/3",
      lambda a, b, p: (11 * b - a) / (6 * b), "(11β-α)/(6β)"),
    R(2, COASTING, ACM, 1, "νc<ν<νm", _c(0.5), "1/2",
      lambda k, p: (3.0 * k - 8.0) / 4.0, "(3k-8)/4",
      lambda a, b, p: (2 * a + 8 * b) / (3 * b), "(2α+8β)/(3β)"),
    R(2, COASTING, ACM, 1, "νm<ν", _p_over_2, "p/2",
      lambda k, p: (2.0 * k + k * p - 8.0) / 4.0, "(2k+kp-8)/4",
      lambda a, b, p: (2 * a + 4) / (1 + b), "(2α+4)/(1+β)"),

    # ---- Table 2, block 2: νa<νm<νc (slow cooling) ----
    R(2, COASTING, AMC, 2, "ν<νa", _c(-2.0), "-2",
      lambda k, p: -2.0, "-2",
      None, "−", note=_T2_AMC_NOTE),
    R(2, COASTING, AMC, 2, "νa<ν<νm", _c(-1.0 / 3.0), "-1/3",
      lambda k, p: -(9.0 - 4.0 * k) / 3.0, "-(9-4k)/3",
      lambda a, b, p: (9 * b - a) / (4 * b), "(9β-α)/(4β)"),
    R(2, COASTING, AMC, 2, "νm<ν<νc", _p_half, "(p-1)/2",
      lambda k, p: (5.0 * k + k * p - 12.0) / 4.0, "(5k+kp-12)/4",
      lambda a, b, p: (2 * a + 6) / (3 + b), "(2α+6)/(3+β)"),
    R(2, COASTING, AMC, 2, "νc<ν", _p_over_2, "p/2",
      lambda k, p: (2.0 * k + k * p - 8.0) / 4.0, "(2k+kp-8)/4",
      lambda a, b, p: (2 * a + 4) / (1 + b), "(2α+4)/(1+β)"),

    # ---- Table 2, block 3: νm<νa<νc (slow cooling) ----
    R(2, COASTING, MAC, 3, "ν<νm", _c(-2.0), "-2",
      lambda k, p: -2.0, "-2",
      None, "−", note=_T2_AMC_NOTE),
    R(2, COASTING, MAC, 3, "νm<ν<νa", _c(-2.5), "-5/2",
      lambda k, p: -(8.0 + k) / 4.0, "-(8+k)/4",
      lambda a, b, p: (10 * a - 8 * b) / b, "(10α-8β)/β"),
    R(2, COASTING, MAC, 3, "νa<ν<νc", _p_half, "(p-1)/2",
      lambda k, p: (5.0 * k + k * p - 12.0) / 4.0, "(5k+kp-12)/4",
      lambda a, b, p: (2 * a + 6) / (3 + b), "(2α+6)/(3+β)"),
    R(2, COASTING, MAC, 3, "νc<ν", _p_over_2, "p/2",
      lambda k, p: (2.0 * k + k * p - 8.0) / 4.0, "(2k+kp-8)/4",
      lambda a, b, p: (2 * a + 4) / (1 + b), "(2α+4)/(1+β)"),

    # ========================= Table 3: BM ========================= #
    # ---- block 1: νa<νc<νm (fast cooling) ----
    R(3, BM, ACM, 1, "ν<νa", _c(-2.0), "-2",
      lambda k, p: -4.0 / (4.0 - k), "-4/(4-k)",
      lambda a, b, p: (4 * a - 2 * b) / a, "(4α-2β)/α"),
    R(3, BM, ACM, 1, "νa<ν<νc", _c(-1.0 / 3.0), "-1/3",
      lambda k, p: -(3.0 * k - 2.0) / (3.0 * (k - 4.0)), "-(3k-2)/[3(k-4)]",
      lambda a, b, p: (4 * a - 2 * b) / (a - 3 * b), "(4α-2β)/(α-3β)"),
    R(3, BM, ACM, 1, "νc<ν<νm", _c(0.5), "1/2",
      lambda k, p: 0.25, "1/4",
      None, "−"),
    R(3, BM, ACM, 1, "νm<ν", _p_over_2, "p/2",
      lambda k, p: (3.0 * p - 2.0) / 4.0, "(3p-2)/4",
      None, "−"),

    # ---- block 2: νa<νm<νc (slow cooling) ----
    R(3, BM, AMC, 2, "ν<νa", _c(-2.0), "-2",
      lambda k, p: -2.0 / (4.0 - k), "-2/(4-k)",
      lambda a, b, p: (4 * a - b) / a, "(4α-β)/α"),
    R(3, BM, AMC, 2, "νa<ν<νm", _c(-1.0 / 3.0), "-1/3",
      lambda k, p: -(k - 2.0) / (k - 4.0), "-(k-2)/(k-4)",
      lambda a, b, p: (4 * a - 6 * b) / (a - 3 * b), "(4α-6β)/(α-3β)"),
    R(3, BM, AMC, 2, "νm<ν<νc", _p_half, "(p-1)/2",
      lambda k, p: (12.0 * p - 12.0 + 5.0 * k - 3.0 * k * p) / (16.0 - 4.0 * k),
      "(12p-12+5k-3kp)/(16-4k)",
      lambda a, b, p: (8 * a - 12 * b) / (2 * a - 3 * b + 1), "(8α-12β)/(2α-3β+1)"),
    R(3, BM, AMC, 2, "νc<ν", _p_over_2, "p/2",
      lambda k, p: (3.0 * p - 2.0) / 4.0, "(3p-2)/4",
      None, "−"),

    # ---- block 3: νm<νa<νc (slow cooling) ----
    R(3, BM, MAC, 3, "ν<νm", _c(-2.0), "-2",
      lambda k, p: -2.0 / (4.0 - k), "-2/(4-k)",
      lambda a, b, p: (4 * a - b) / a, "(4α-β)/α"),
    R(3, BM, MAC, 3, "νm<ν<νa", _c(-2.5), "-5/2",
      lambda k, p: -(20.0 - 3.0 * k) / (4.0 * (4.0 - k)), "-(20-3k)/[4(4-k)]",
      lambda a, b, p: (40 * a - 20 * b) / (10 * a - 3 * b), "(40α-20β)/(10α-3β)"),
    R(3, BM, MAC, 3, "νa<ν<νc", _p_half, "(p-1)/2",
      lambda k, p: (12.0 * p - 12.0 + 5.0 * k - 3.0 * k * p) / (16.0 - 4.0 * k),
      "(12p-12+5k-3kp)/(16-4k)",
      lambda a, b, p: (8 * a - 12 * b) / (2 * a - 3 * b + 1), "(8α-12β)/(2α-3β+1)"),
    R(3, BM, MAC, 3, "νc<ν", _p_over_2, "p/2",
      lambda k, p: (3.0 * p - 2.0) / 4.0, "(3p-2)/4",
      None, "−"),

    # ========================= Table 4: ST ========================= #
    # ---- block 1: νa<νm<νc ----
    R(4, ST, AMC, 1, "ν<νa", _c(-2.0), "-2",
      lambda k, p: -2.0 * (k - 1.0) / (5.0 - k), "-2(k-1)/(5-k)",
      lambda a, b, p: (5 * a + b) / (a + b), "(5α+β)/(α+β)"),
    R(4, ST, AMC, 1, "νa<ν<νm", _c(-1.0 / 3.0), "-1/3",
      lambda k, p: -(24.0 - 10.0 * k) / (3.0 * (5.0 - k)), "-(24-10k)/[3(5-k)]",
      lambda a, b, p: (5 * a - 24 * b) / (a - 10 * b), "(5α-24β)/(α-10β)"),
    R(4, ST, AMC, 1, "νm<ν<νc", _p_half, "(p-1)/2",
      lambda k, p: (15.0 * p - 4.0 * k * p + 8.0 * k - 21.0) / (2.0 * (5.0 - k)),
      "(15p-4kp+8k-21)/[2(5-k)]",
      lambda a, b, p: (5 * a - 15 * b + 3) / (a - 4 * b + 2), "(5α-15β+3)/(α-4β+2)"),
    R(4, ST, AMC, 1, "νc<ν", _p_over_2, "p/2",
      lambda k, p: (15.0 * p - 4.0 * k * p + 6.0 * k - 20.0) / (2.0 * (5.0 - k)),
      "(15p-4kp+6k-20)/[2(5-k)]",
      lambda a, b, p: (5 * a - 15 * b + 10) / (a - 4 * b + 3), "(5α-15β+10)/(α-4β+3)"),

    # ---- block 2: νm<νa<νc ----
    R(4, ST, MAC, 2, "ν<νm", _c(-2.0), "-2",
      lambda k, p: -2.0 * (k - 1.0) / (5.0 - k), "-2(k-1)/(5-k)",
      lambda a, b, p: (5 * a + b) / (a + b), "(5α+β)/(α+β)"),
    R(4, ST, MAC, 2, "νm<ν<νa", _c(-2.5), "-5/2",
      lambda k, p: -11.0 / (2.0 * (5.0 - k)), "-11/[2(5-k)]",
      lambda a, b, p: (25 * a - 11 * b) / (5 * a), "(25α-11β)/(5α)"),
    R(4, ST, MAC, 2, "νa<ν<νc", _p_half, "(p-1)/2",
      lambda k, p: (15.0 * p - 4.0 * k * p + 8.0 * k - 21.0) / (2.0 * (5.0 - k)),
      "(15p-4kp+8k-21)/[2(5-k)]",
      lambda a, b, p: (5 * a - 15 * b + 3) / (a - 4 * b + 2), "(5α-15β+3)/(α-4β+2)"),
    R(4, ST, MAC, 2, "νc<ν", _p_over_2, "p/2",
      lambda k, p: (15.0 * p - 4.0 * k * p + 6.0 * k - 20.0) / (2.0 * (5.0 - k)),
      "(15p-4kp+6k-20)/[2(5-k)]",
      lambda a, b, p: (5 * a - 15 * b + 10) / (a - 4 * b + 3), "(5α-15β+10)/(α-4β+3)"),
]

# ====================== Table 5: deep-Newtonian =========================== #
_DN_NOTE = ("Table 5 按 νm<νa<νc 推导；对 νa<νm<νc，其 νm<ν<νc 段的 α(k) "
            "与 νa<ν<νc 段相同（与 Table 4 同样的简并），故公式通用")
ROWS += [
    # ---- block 1: 2 < p < 3 ----
    R(5, DN, MAC, 1, "ν<νm", _c(-2.0), "-2",
      lambda k, p: -4.0 / (5.0 - k), "-4/(5-k)",
      lambda a, b, p: (5 * a - 2 * b) / a, "(5α-2β)/α",
      prange=(2.0, 3.0), note=_DN_NOTE),
    R(5, DN, MAC, 1, "νm<ν<νa", _c(-2.5), "-5/2",
      lambda k, p: -11.0 / (2.0 * (5.0 - k)), "-11/[2(5-k)]",
      lambda a, b, p: (25 * a - 11 * b) / (5 * a), "(25α-11β)/(5α)",
      prange=(2.0, 3.0), note=_DN_NOTE),
    R(5, DN, MAC, 1, "νa<ν<νc", _p_half, "(p-1)/2",
      lambda k, p: 3.0 * (1.0 + p) / (2.0 * (5.0 - k)), "3(1+p)/[2(5-k)]",
      lambda a, b, p: (5 * a - 3 * b - 3) / a, "(5α-3β-3)/α",
      prange=(2.0, 3.0), note=_DN_NOTE),
    R(5, DN, MAC, 1, "νc<ν", _p_over_2, "p/2",
      lambda k, p: (3.0 * p + 4.0 - 2.0 * k) / (2.0 * (5.0 - k)), "(3p+4-2k)/[2(5-k)]",
      lambda a, b, p: (5 * a - 3 * b - 2) / (a - 1), "(5α-3β-2)/(α-1)",
      prange=(2.0, 3.0), note=_DN_NOTE),

    # ---- block 2: p ≳ 3 ----
    R(5, DN, MAC, 2, "ν<νm", _c(-2.0), "-2",
      lambda k, p: -(k * p - 3.0 * p - 3.0 * k + 13.0) / (5.0 - k),
      "-(kp-3p-3k+13)/(5-k)",
      lambda a, b, p: (10 * a - 13 * b + 3 * p * b) / (2 * a - 3 * b + p * b),
      "(10α-13β+3pβ)/(2α-3β+pβ)",
      prange=(3.0, math.inf), note=_DN_NOTE),
    R(5, DN, MAC, 2, "νm<ν<νa", _c(-2.5), "-5/2",
      lambda k, p: -(2.0 * k * p - 6.0 * p - 6.0 * k + 29.0) / (2.0 * (5.0 - k)),
      "-(2kp-6p-6k+29)/[2(5-k)]",
      lambda a, b, p: (25 * a - 29 * b + 6 * p * b) / (5 * a - 6 * b + 2 * p * b),
      "(25α-29β+6pβ)/(5α-6β+2pβ)",
      prange=(3.0, math.inf), note=_DN_NOTE),
    R(5, DN, MAC, 2, "νa<ν<νc", _p_half, "(p-1)/2",
      lambda k, p: (9.0 * p + 6.0 * k - 2.0 * k * p - 15.0) / (2.0 * (5.0 - k)),
      "(9p+6k-2kp-15)/[2(5-k)]",
      lambda a, b, p: (5 * a - 9 * b + 3) / (a - 2 * b + 2), "(5α-9β+3)/(α-2β+2)",
      prange=(3.0, math.inf), note=_DN_NOTE),
    R(5, DN, MAC, 2, "νc<ν", _p_over_2, "p/2",
      lambda k, p: (9.0 * p + 4.0 * k - 2.0 * k * p - 14.0) / (2.0 * (5.0 - k)),
      "(9p+4k-2kp-14)/[2(5-k)]",
      lambda a, b, p: (5 * a - 9 * b + 7) / (a - 2 * b + 2), "(5α-9β+7)/(α-2β+2)",
      prange=(3.0, math.inf), note=_DN_NOTE),
]

# 按表号/排序索引
TABLE: Dict[int, Dict[str, List[Row]]] = {}
for _r in ROWS:
    TABLE.setdefault(_r.table, {}).setdefault(_r.ordering, []).append(_r)


# --------------------------------------------------------------------------- #
# 2. 选择与计算
# --------------------------------------------------------------------------- #
def assumed_beta(outflow: str, p: float) -> float:
    """按需求规则给出 β：jetted → -2；非相对论 → (p-1)/2。"""
    if outflow == "jetted":
        return -2.0
    if outflow == "nonrel":
        return (p - 1.0) / 2.0
    raise ValueError(tr("err_outflow", v=outflow))


def phase_table(phase: str) -> int:
    return {COASTING: 2, BM: 3, ST: 4, DN: 5}[phase]


_PHASE_ALIAS: Dict[str, str] = {
    "coasting": COASTING, "bm": BM, "st": ST, "dn": DN,
}


def norm_phase(phase: str) -> str:
    """把 CLI/交互里的小写写法归一成内部常量。"""
    v = _PHASE_ALIAS.get(str(phase).lower())
    if v is None:
        raise ValueError(tr("err_phase", v=phase))
    return v


def select_row(table: int, ordering: Optional[str], beta: float, p: float) -> Row:
    """在指定表中选出 β 与之匹配的唯一一行。"""
    if table == 5:
        order_key = MAC          # Table 5 仅按 νm<νa<νc 推导
    else:
        order_key = ordering or AMC
    rows = TABLE.get(table, {}).get(order_key, [])
    cand = [r for r in rows if r.in_p_range(p) and r.match_beta(beta, p)]

    if len(cand) == 1:
        return cand[0]
    if not cand:
        raise ValueError(tr("err_row_not_found", tbl=TABLE_NAME[table],
                            ord=ORDER_LABEL_PLAIN.get(order_key, order_key),
                            b=f"{beta:g}", p=f"{p:g}"))
    raise ValueError(tr("err_row_ambiguous", tbl=TABLE_NAME[table],
                        b=f"{beta:g}", n=len(cand),
                        segs="; ".join(r.segment for r in cand)))


def _safe_k(row: Row, alpha: float, beta: float, p: float) -> Tuple[Optional[float], Optional[str]]:
    """计算 k_CR，处理表中为 '−' 及奇异（分母为零/非有限值）的情况。"""
    if row.kcr is None:
        return None, None
    try:
        k = float(row.kcr(alpha, beta, p))
    except ZeroDivisionError:
        return None, tr("msg_zero_den")
    if not math.isfinite(k):
        return None, tr("msg_not_finite")
    return k, None


def _num_partial(f: Callable[[float], float], x: float, h: float = 1e-5) -> Optional[float]:
    """中心差分求导；奇异时返回 None。"""
    try:
        return (f(x + h) - f(x - h)) / (2.0 * h)
    except (ZeroDivisionError, ValueError):
        return None


def _propagate(row: Row, alpha: float, p: float, outflow: str,
               aerr: Optional[float], perr: Optional[float]) -> Optional[float]:
    """误差传递（数值偏导）：δk = sqrt[(∂k/∂α δα)^2 + (∂k/∂p δp)^2]。"""
    if row.kcr is None or (not aerr and not perr):
        return None

    def k_of(a: float, pp: float) -> float:
        return row.kcr(a, assumed_beta(outflow, pp), pp)

    parts = []
    if aerr:
        d = _num_partial(lambda a: k_of(a, p), alpha)
        if d is not None:
            parts.append((d * aerr) ** 2)
    if perr:
        d = _num_partial(lambda pp: k_of(alpha, pp), p)
        if d is not None:
            parts.append((d * perr) ** 2)
    if not parts:
        return None
    return math.sqrt(sum(parts))


@dataclass
class Result:
    alpha: float
    p: float
    outflow: str
    cooling: Optional[str]
    ordering: Optional[str]
    phase: str
    beta: float
    row: Row
    k: Optional[float]
    k_err: Optional[float]
    alpha_model: float
    warnings: List[str] = field(default_factory=list)
    path: str = ""
    message: str = ""
    alpha_err: Optional[float] = None
    p_err: Optional[float] = None


def evaluate(alpha: float, p: float, outflow: str, phase: str,
             cooling: Optional[str] = None, ordering: Optional[str] = None,
             alpha_err: Optional[float] = None,
             p_err: Optional[float] = None) -> Result:
    """按规则选行并计算 k_CR。"""
    phase = norm_phase(phase)
    if ordering is not None:
        ordering = str(ordering).upper()
        if ordering not in (ACM, AMC, MAC):
            raise ValueError(tr("err_ordering", v=ordering))
    beta = assumed_beta(outflow, p)

    # 选表 + 选 block
    if outflow == "jetted":
        order_eff = ACM if cooling == "fast" else (ordering or AMC)
    else:
        order_eff = ordering or AMC

    table = phase_table(phase)
    row = select_row(table, order_eff, beta, p)

    k, err_msg = _safe_k(row, alpha, beta, p)
    if row.kcr is None:
        message = tr("msg_kcr_none", lab=row.alpha_label)
    else:
        message = err_msg or ""
    k_err = _propagate(row, alpha, p, outflow, alpha_err, p_err) if k is not None else None
    alpha_model = row.alpha_of_k(k if k is not None else 0.0, p)

    warnings: List[str] = []
    if err_msg:
        warnings.append(err_msg)
    if p <= 2.0:
        warnings.append(tr("w_p_le2", p=f"{p:g}"))
    elif p > 4.0:
        warnings.append(tr("w_p_gt4", p=f"{p:g}"))
    if phase == DN and abs(p - 3.0) < 0.15:
        warnings.append(tr("w_p_near3", p=f"{p:g}"))
    if k is not None:
        if k < 0.0:
            warnings.append(tr("w_k_neg", k=f"{k:.3f}"))
        elif k > 3.0:
            warnings.append(tr("w_k_gt3", k=f"{k:.3f}"))
        else:
            warnings.append(tr("w_k_ok", k=f"{k:.3f}"))
        # α 是否落在 k∈[0,3] 对应的区间内（判断外插）
        a_lo, a_hi = sorted((row.alpha_of_k(0.0, p), row.alpha_of_k(3.0, p)))
        if not (a_lo - 1e-9 <= alpha <= a_hi + 1e-9):
            warnings.append(tr("w_alpha_extrap", a=f"{alpha:g}",
                               lo=f"{a_lo:.3f}", hi=f"{a_hi:.3f}"))
    elif row.kcr is None:
        warnings.append(tr("w_kcr_none", lab=row.alpha_label))

    # 动力学阶段的先验提示
    if outflow == "jetted" and phase == COASTING:
        warnings.append(tr("w_jetted_coasting"))
    if outflow == "nonrel" and phase != COASTING:
        for r in candidate_rows(outflow, p):
            if r.table == 2 and r.kcr is not None:
                kc, _ = _safe_k(r, alpha, beta, p)
                if kc is not None and 0.0 <= kc <= 3.0:
                    warnings.append(tr("w_nonrel_alt", k=f"{kc:.2f}"))
                    break
    if phase == DN and order_eff == AMC:
        warnings.append(tr("w_dn_amc"))

    path = _path_text(outflow, cooling, order_eff, phase, row)
    return Result(alpha, p, outflow, cooling, order_eff, phase, beta, row, k, k_err,
                  alpha_model, warnings, path, message, alpha_err, p_err)


def _path_text(outflow: str, cooling: Optional[str], ordering: str,
               phase: str, row: Row) -> str:
    if outflow == "jetted":
        steps = [tr("path_jetted"), "β=-2"]
        steps.append(tr("path_fast") if cooling == "fast" else tr("path_slow"))
        steps.append(ORDER_LABEL_PLAIN[ordering])
        steps.append(phase_label(phase))
    else:
        steps = [tr("path_nonrel"), "β=(p-1)/2", ORDER_LABEL_PLAIN[ordering],
                 phase_label(phase)]
    return " → ".join(steps)


# --------------------------------------------------------------------------- #
# 3. 诊断：同一 β 下其它动力学阶段的 k_CR
# --------------------------------------------------------------------------- #
def candidate_rows(outflow: str, p: float) -> List[Row]:
    """同一 β 下、适用该外流类型的其它行（用于对比）。

    jetted  : Table 2 (coasting) + Table 3 (BM)
    非相对论 : Table 2 (coasting) + Table 4 (ST) + Table 5 (DN)
               （Table 3 的 BM 只适用于相对论喷流，故不列出）
    """
    beta = assumed_beta(outflow, p)
    allowed = {2, 3} if outflow == "jetted" else {2, 4, 5}
    out = []
    for row in ROWS:
        if row.table not in allowed or not row.in_p_range(p):
            continue
        if row.match_beta(beta, p):
            out.append(row)
    return out


# --------------------------------------------------------------------------- #
# 4. 输出
# --------------------------------------------------------------------------- #
BAR = "─" * 78
DBL = "═" * 78


def _f(x: Optional[float], nd: int = 3) -> str:
    return "—" if x is None else f"{x:.{nd}f}"


def report(res: Result, show_candidates: bool = True) -> str:
    L: List[str] = []
    L.append(DBL)
    L.append(tr("rep_banner"))
    L.append(DBL)
    a_txt = f"{res.alpha:g}" + (f" +- {res.alpha_err:g}" if res.alpha_err else "")
    p_txt = f"{res.p:g}" + (f" +- {res.p_err:g}" if res.p_err else "")
    L.append(tr("rep_input", a=a_txt, b=p_txt))
    L.append(tr("rep_path", p=res.path))
    L.append(tr("rep_selected", w=res.row.where()))
    L.append(tr("rep_beta", b=f"{res.beta:g}", lab=res.row.beta_label))
    L.append(BAR)
    L.append(tr("rep_alphak", f=res.row.alpha_label))
    L.append(tr("rep_kcrformula", f=res.row.kcr_label))
    L.append(BAR)
    if res.k is not None:
        err = f" +- {res.k_err:.3f}" if res.k_err is not None else ""
        L.append(tr("rep_k", k=f"{res.k:.4f}", err=err))
        resid = res.alpha_model - res.alpha
        L.append(tr("rep_selfcheck", m=f"{res.alpha_model:.4f}",
                    a=f"{res.alpha:.4f}", r=f"{resid:+.1e}"))
    else:
        L.append(tr("rep_undet"))
        L.append(f"        {res.message}")
        if res.row.kcr is None:
            a_pred = res.row.alpha_of_k(0.0, res.p)   # 该行 α 与 k 无关
            key = "tip_close" if abs(res.alpha - a_pred) < 0.3 else "tip_off"
            L.append(tr(key, p=f"{a_pred:g}", a=f"{res.alpha:g}"))
        else:
            a0 = res.row.alpha_of_k(0.0, res.p)
            a3 = res.row.alpha_of_k(3.0, res.p)
            L.append(tr("tip_singular", a0=f"{a0:.2f}", a3=f"{a3:.2f}",
                        a=f"{res.alpha:g}"))
    if res.warnings:
        L.append(BAR)
        L.append(tr("rep_checks"))
        for w in res.warnings:
            L.append(f"  • {w}")

    if show_candidates:
        L.append(BAR)
        L.append(tr("rep_candidates", b=f"{res.beta:g}"))
        for row in candidate_rows(res.outflow, res.p):
            k, _ = _safe_k(row, res.alpha, res.beta, res.p)
            mark = tr("rep_mark") if row == res.row else ""
            L.append(f"  {TABLE_NAME[row.table]:<30s} block{row.block} "
                     f"{ORDER_LABEL_PLAIN[row.ordering]:<12s} {row.segment:<8s} "
                     f"k = {_f(k):>9s}   [{row.kcr_label}]{mark}")
        L.append(tr("rep_foot_jetted") if res.outflow == "jetted"
                 else tr("rep_foot_nonrel"))
    L.append(DBL)
    return "\n".join(L)


def result_to_dict(res: Result) -> dict:
    cands = []
    for row in candidate_rows(res.outflow, res.p):
        k, msg = _safe_k(row, res.alpha, res.beta, res.p)
        cands.append({
            "table": row.table, "phase": row.phase, "ordering": row.ordering,
            "block": row.block, "segment": row.segment,
            "beta_label": row.beta_label, "alpha_k": row.alpha_label,
            "kcr_formula": row.kcr_label,
            "k": k, "error": msg, "selected": row == res.row,
        })
    return {
        "input": {"alpha": res.alpha, "p": res.p,
                  "outflow": res.outflow, "cooling": res.cooling,
                  "ordering": res.ordering, "phase": res.phase},
        "beta": res.beta,
        "selected": {"table": res.row.table, "phase": res.row.phase,
                     "ordering": res.row.ordering, "block": res.row.block,
                     "segment": res.row.segment,
                     "alpha_k": res.row.alpha_label,
                     "kcr_formula": res.row.kcr_label,
                     "where": res.row.where()},
        "k_CR": res.k,
        "k_CR_err": res.k_err,
        "alpha_model": res.alpha_model,
        "warnings": res.warnings,
        "candidates": cands,
    }


# --------------------------------------------------------------------------- #
# 5. 自检：用 α(k) 反解 k_CR，验证 Table 2-5 每行公式
# --------------------------------------------------------------------------- #
def selftest(verbose: bool = True) -> bool:
    k_grid = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    ok = True
    n_pass = n_fail = n_skip = 0

    for row in ROWS:
        # β 为 p/2 或 (p-1)/2 时，挑使 β 命中表格数值的 p
        if row.beta_label in ("p/2", "(p-1)/2"):
            bpn = row.beta_fn(1.0)
            need = 1.0 if abs(bpn - 0.0) < 1e-12 else 0.5
            p_list = [pp for pp in [2.0 + 2 * need, 2.5, 3.5, 4.0]]
        else:
            p_list = [2.1, 2.5, 3.0, 3.5]

        for p in p_list:
            if not row.in_p_range(p):
                continue
            beta = row.beta(p)
            if row.kcr is None:
                # 表中 "−"：α 必须与 k 无关
                alphas = [row.alpha_of_k(k, p) for k in k_grid]
                if max(alphas) - min(alphas) > 1e-12:
                    print(tr("st_fail_none", w=row.where()))
                    ok = False
                    n_fail += 1
                else:
                    n_pass += 1
                continue
            for k in k_grid:
                alpha = row.alpha_of_k(k, p)
                k_back, _ = _safe_k(row, alpha, beta, p)
                if k_back is None:
                    n_skip += 1
                    continue
                if abs(k_back - k) > 1e-7 * max(1.0, abs(k)):
                    print(tr("st_fail_inv", w=row.where(), p=f"{p:g}", k=f"{k:g}",
                             a=f"{alpha:.10g}", kb=f"{k_back:.10g}"))
                    ok = False
                    n_fail += 1
                else:
                    n_pass += 1

    if verbose:
        print(tr("st_summary", p=n_pass, f=n_fail, s=n_skip))
        print(tr("st_ok") if ok else tr("st_bad"))
    return ok


# --------------------------------------------------------------------------- #
# 6. 用论文 Table 6 的 26 个源复核（label → 表/block/行）
# --------------------------------------------------------------------------- #
# (name, β_CR, p, α_CR, k_CR_table, label)
PAPER_TABLE6: List[Tuple[str, float, float, float, float, str]] = [
    ("Swift J1644+57",          -2.00, 2.28, -1.16, 2.27, "a1"),
    ("AT2022cmc",               -2.00, 2.15, -0.72, 1.23, "a1"),
    ("ASASSN-15oi",              1.15, 3.31,  3.70, 3.22, "b1"),
    ("AT2024tvd",                1.43, 3.85,  4.37, 3.33, "b1"),
    ("AT2020vwl",                1.47, 3.95, -0.61, 1.07, "b1"),
    ("ASASSN-14li",              1.30, 3.60,  1.68, 2.18, "b1"),
    ("AT2020wjw",                1.47, 3.94, -0.04, 1.32, "b2"),
    ("ASASSN-14ae",              0.51, 2.03, -5.13, -1.21, "b1"),
    ("AT2018hyz",                0.50, 2.01, -4.20, -0.69, "b2"),
    ("AT2019ahk",                1.26, 3.52, -0.34, 1.25, "b2"),
    ("AT2019avd",                1.07, 3.14,  0.00, 1.48, "b2"),
    ("AT2019azh",                1.50, 3.99,  0.23, 1.43, "b1"),
    ("AT2019dsg",                0.86, 2.72,  0.01, 1.56, "b2"),
    ("AT2019ehz",                1.46, 3.92,  3.81, 3.05, "b2"),
    ("AT2019eve",                1.42, 3.85,  0.26, 1.47, "b1"),
    ("AT2019teq",                0.85, 2.71,  0.18, 1.65, "b2"),
    ("AT2020opy",                1.20, 3.41,  0.01, 1.43, "b2"),
    ("AT2020vdq",                1.05, 3.10,  0.03, 1.50, "b2"),
    ("CNSS J0019+00",            1.02, 3.04,  2.55, 2.76, "b2"),
    ("eRASSt J011431-593654",    1.30, 3.59,  0.09, 1.44, "b1"),
    ("eRASSt J210858-562832",    1.15, 3.31,  0.08, 1.48, "b2"),
    ("IGR J12580+0134",          1.49, 3.99,  0.01, 1.34, "b2"),
    ("iPTF 16fnl",               1.13, 3.26,  0.09, 1.49, "b2"),
    ("PS16dtm",                  1.46, 3.91,  0.21, 1.44, "b1"),
    ("WTP14adeqka",              1.18, 3.36, -1.98, 0.49, "b1"),
    ("XMMSL1 J0740-85",          1.47, 3.95,  0.31, 1.48, "b1"),
]
# 论文 Table 6 脚注：a=BM 减速（a1: νa<νm<νc, a2: νm<νa<νc）
#                b=coasting（b1: νa<νm<νc, b2: νm<νa<νc）
_LABEL_TO_ROW: Dict[str, Tuple[int, str, str]] = {
    "a1": (3, AMC, "ν<νa"),      # BM, β=-2, 最低频段
    "a2": (3, MAC, "ν<νm"),
    "b1": (2, AMC, "νm<ν<νc"),   # coasting, β=(p-1)/2
    "b2": (2, MAC, "νa<ν<νc"),
}


def verify_paper(verbose: bool = True) -> bool:
    ok = True
    rows_out: List[str] = []
    n_pass = 0
    for name, beta, p, alpha, k_tab, label in PAPER_TABLE6:
        table, ordering, seg = _LABEL_TO_ROW[label]
        row = next(r for r in TABLE[table][ordering] if r.segment == seg)
        k_pred, _ = _safe_k(row, alpha, beta, p)
        good = k_pred is not None and abs(k_pred - k_tab) <= 0.02
        n_pass += int(good)
        ok &= good
        d = _f(None if k_pred is None else k_pred - k_tab, 3)
        rows_out.append(
            f"  {name:<24s} {label}  β={beta:>5.2f} α={alpha:>6.2f}  "
            + tr("vp_cols", kt=f"{k_tab:>6.2f}", kc=f"{_f(k_pred):>6s}", d=f"{d:>7s}")
            + f"  {'OK' if good else 'MISMATCH'}")
    if verbose:
        print(tr("vp_header"))
        print(tr("vp_legend"))
        for line in rows_out:
            print(line)
        print(tr("vp_summary", n=n_pass, N=len(PAPER_TABLE6)))
    return ok


# --------------------------------------------------------------------------- #
# 7. 打印 Table 2-5
# --------------------------------------------------------------------------- #
def list_tables() -> None:
    for tno in (2, 3, 4, 5):
        print(DBL)
        print(f" {TABLE_NAME[tno]}")
        print(DBL)
        for ordering, rows in TABLE[tno].items():
            for i, row in enumerate(rows):
                head = (f"  block{row.block}  {ORDER_LABEL_PLAIN[ordering]:<12s} "
                        f"{row.segment:<10s} β={row.beta_label:<8s}") if i == 0 or \
                       rows[i - 1].block != row.block else \
                       (f"  {'':<6s}  {ORDER_LABEL_PLAIN[ordering]:<12s} "
                        f"{row.segment:<10s} β={row.beta_label:<8s}")
                print(f"{head} α(k)={row.alpha_label:<26s} k_CR={row.kcr_label}")
            print()
    return


# --------------------------------------------------------------------------- #
# 8. 交互式
# --------------------------------------------------------------------------- #
_UNICODE_MAP = str.maketrans({
    "−": "-", "–": "-", "—": "-", "﹣": "-",     # 各种减号
    "＋": "+", "＝": "=", "．": ".", "。": ".",
})


def parse_value(text: str) -> Tuple[float, Optional[float]]:
    """解析 "值" 或 "值+误差"：误差分隔符支持 +- 与 -+（也兼容 +/-）。"""
    s = text.translate(_UNICODE_MAP).strip().replace(" ", "")
    for sep in ("+-", "-+", "+/-", "±"):
        if sep in s:
            a, e = s.split(sep, 1)
            return float(a), abs(float(e))
    return float(s), None


def ask_value(prompt: str, default: Optional[float] = None) -> Tuple[float, Optional[float]]:
    while True:
        try:
            s = input(prompt).strip()
        except EOFError:
            if default is None:
                raise
            return default, None
        if not s and default is not None:
            return default, None
        try:
            return parse_value(s)
        except ValueError:
            print(tr("it_bad_number"))


def ask_choice(prompt: str, options: Sequence[Tuple[str, str]]) -> str:
    keys = [k for k, _ in options]
    menu = prompt + "   " + "  ".join(
        f"[{i}] {lab}" for i, (k, lab) in enumerate(options, 1))
    while True:
        try:
            s = input(menu + " > ").strip().lower()
        except EOFError:
            return keys[0]
        if s.isdigit() and 1 <= int(s) <= len(options):
            return options[int(s) - 1][0]
        for k in keys:
            if s == k.lower():
                return k
        print(tr("it_bad_choice"))


def interactive() -> int:
    print(tr("it_title"))
    print(tr("it_err_fmt"))
    alpha, aerr = ask_value(tr("it_alpha"))
    p, perr = ask_value(tr("it_p"))

    outflow = ask_choice(tr("it_ask_outflow"), [
        ("jetted", tr("it_opt_jetted")),
        ("nonrel", tr("it_opt_nonrel")),
    ])

    cooling = ordering = None
    phases: List[str] = []

    if outflow == "jetted":
        cooling = ask_choice(tr("it_ask_cooling"), [
            ("fast", tr("it_opt_fast")),
            ("slow", tr("it_opt_slow")),
        ])
        if cooling == "fast":
            ordering = ACM
            phases = [ask_choice(tr("it_ask_phase", n=3), [
                ("coasting", tr("it_opt_coasting")),
                ("bm", tr("it_opt_bm")),
            ])]
        else:
            ordering = ask_choice(tr("it_ask_order", n=3), [
                ("amc", tr("it_opt_amc")),
                ("mac", tr("it_opt_mac")),
            ])
            phases = [ask_choice(tr("it_ask_phase", n=4), [
                ("coasting", tr("it_opt_coasting")),
                ("bm", tr("it_opt_bm")),
            ])]
    else:
        ordering = ask_choice(tr("it_ask_order", n=2), [
            ("amc", tr("it_opt_amc")),
            ("mac", tr("it_opt_mac")),
        ])
        if alpha < 1.0:
            print(tr("it_alpha_lt1", a=f"{alpha:g}"))
            phases = [COASTING]
        else:
            phase = ask_choice(tr("it_ask_decel"), [
                ("st", tr("it_opt_st")),
                ("dn", tr("it_opt_dn")),
                ("both", tr("it_opt_both")),
            ])
            phases = [ST, DN] if phase == "both" else [phase]

    outs = []
    for ph in phases:
        res = evaluate(alpha, p, outflow, ph, cooling, ordering, aerr, perr)
        print()
        print(report(res))
        outs.append(result_to_dict(res))
    return 0


# --------------------------------------------------------------------------- #
# 9. CLI
# --------------------------------------------------------------------------- #
def printf_all(args) -> int:
    alpha, aerr = parse_value(args.alpha)
    p, perr = parse_value(args.p)

    outflow = "jetted" if args.jetted else ("nonrel" if args.nonrel else None)
    if outflow is None:
        raise SystemExit(tr("cli_need_outflow"))

    cooling = args.cooling
    ordering = args.order.upper() if args.order else None
    phase_arg = args.phase

    if outflow == "jetted":
        if cooling is None:
            raise SystemExit(tr("cli_need_cooling"))
        if phase_arg is None:
            raise SystemExit(tr("cli_need_phase"))
        if cooling == "fast" and ordering not in (None, ACM):
            print(tr("cli_note_fast"), file=sys.stderr)
        if cooling == "slow" and ordering is None:
            raise SystemExit(tr("cli_need_order_slow"))
        phases = [phase_arg]
    else:
        if ordering is None:
            raise SystemExit(tr("cli_need_order_nonrel"))
        if phase_arg is None:
            phase_arg = "coasting" if alpha < 1.0 else None
            if phase_arg is None:
                raise SystemExit(tr("cli_need_phase_alpha1"))
        phases = [ST, DN] if phase_arg == "both" else [phase_arg]

    out = []
    for ph in phases:
        res = evaluate(alpha, p, outflow, ph, cooling, ordering, aerr, perr)
        d = result_to_dict(res)
        out.append(d)
        if not args.json:
            print(report(res))
            print()
    if args.json:
        print(json.dumps(out if len(out) > 1 else out[0], ensure_ascii=False, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=tr("ap_desc"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=tr("ap_epilog"))
    ap.add_argument("-a", "--alpha", type=str, help=tr("ap_alpha"))
    ap.add_argument("-p", "--p", dest="p", type=str, help=tr("ap_p"))
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--jetted", action="store_true", help=tr("ap_jetted"))
    g.add_argument("--nonrel", action="store_true", help=tr("ap_nonrel"))
    ap.add_argument("--cooling", choices=["fast", "slow"], help=tr("ap_cooling"))
    ap.add_argument("--order", choices=["acm", "amc", "mac"], help=tr("ap_order"))
    ap.add_argument("--phase", choices=["coasting", "bm", "st", "dn", "both"],
                    help=tr("ap_phase"))
    ap.add_argument("--json", action="store_true", help=tr("ap_json"))
    ap.add_argument("--selftest", action="store_true", help=tr("ap_selftest"))
    ap.add_argument("--verify-paper", action="store_true", help=tr("ap_verify"))
    ap.add_argument("--list-tables", action="store_true", help=tr("ap_list"))
    ap.add_argument("-e", "--english", action="store_true", help=tr("ap_english"))
    return ap


_VALUE_OPTS = {"-a", "--alpha", "-p", "--p"}


def _normalize_argv(argv: Sequence[str]) -> List[str]:
    """让 `-a -1.16+-0.02` 这类以负号开头的取值也能被 argparse 接受。"""
    out: List[str] = []
    i = 0
    while i < len(argv):
        tok = argv[i]
        if (tok in _VALUE_OPTS and i + 1 < len(argv)
                and argv[i + 1].startswith("-") and argv[i + 1] != "-"):
            out.append(f"{tok}={argv[i + 1]}")
            i += 2
            continue
        out.append(tok)
        i += 1
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    raw = list(sys.argv[1:] if argv is None else argv)
    # 语言需在 build_parser 之前确定，否则 --help 的文案不会切换
    set_language("en" if ("-e" in raw or "--english" in raw) else "zh")
    args = build_parser().parse_args(_normalize_argv(raw))
    if args.english:
        set_language("en")

    if args.list_tables:
        list_tables()
        return 0
    if args.selftest:
        return 0 if selftest() else 1
    if args.verify_paper:
        return 0 if verify_paper() else 1

    if args.alpha is None and args.p is None and not (args.jetted or args.nonrel):
        return interactive()

    if args.alpha is None or args.p is None:
        raise SystemExit(tr("cli_need_alpha_p"))
    return printf_all(args)


if __name__ == "__main__":
    sys.exit(main())
