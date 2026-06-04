"""
q ≡ 3 (mod 4) square-field verification / witness search.

This script checks odd prime-power field orders Q < 553736 with Q ≡ 3 (mod 4),
matching the theoretical bound used in the paper. It depends on
pythagorean_fields_common.py and writes progress/results to q3_results.txt.

Default mode:
    b3_then_auto_line

This first tries the very fast center-one line u=t, v=3t. If that line is
invalid or fails, it tries other valid center-one lines u=t, v=b t.

Important interpretation:
    - A found witness proves F_Q is non-Parker.
    - A failure in line-search mode only means no line witness was found.
    - To certify Parker behavior for a failure, rerun those Q with MODE="full_uv".

For Q ≡ 3 (mod 4), center-zero cannot give a distinct square magic square
because -1 is not a square. Therefore a complete center-one full_uv failure
certifies Parker for that Q.

Run:
    python pythagorean_fields_q3.py

Optional environment overrides, for convenience:
    BOUND=1000 MODE=b3_then_auto_line python pythagorean_fields_q3.py
    ORDER_MODE=manual MANUAL_ORDERS=3,7,11,19,23,27,31,43,47,67,243 MODE=full_uv python pythagorean_fields_q3.py
"""

from __future__ import annotations

import gc
import os
from time import perf_counter

from pythagorean_fields_common import (
    FiniteField,
    prime_power_data,
    prime_power_orders,
    center_one_entries,
    find_center_one_full_witness,
    make_witness,
    witness_to_text,
)


# ============================================================
# SETTINGS
# ============================================================

BOUND = int(os.environ.get("BOUND", "553736"))
EXPONENT_N = 2

# Options:
#   "Q3mod4"     : all prime powers Q < BOUND with Q ≡ 3 mod 4
#   "basep3mod4" : all prime powers Q=p^k < BOUND with p ≡ 3 mod 4
#   "manual"     : use MANUAL_ORDERS below
ORDER_MODE = os.environ.get("ORDER_MODE", "Q3mod4")

# Used only when ORDER_MODE="manual".
_default_manual = "3,7,11,19,23,27,31,43,47,67,243"
MANUAL_ORDERS = [int(x) for x in os.environ.get("MANUAL_ORDERS", _default_manual).split(",") if x.strip()]

# Options:
#   "b3_line"            : only the fast line u=t, v=3t
#   "auto_line"          : try valid b-lines u=t, v=b t
#   "b3_then_auto_line"  : try b=3 first, then auto_line
#   "full_uv"            : complete center-one search over all U,V
MODE = os.environ.get("MODE", "b3_then_auto_line")

RESULTS_TXT = os.environ.get("RESULTS_TXT", "q3_results.txt")
GC_EVERY = int(os.environ.get("GC_EVERY", "20"))


# ============================================================
# ORDER SELECTION
# ============================================================


def relevant_orders():
    if ORDER_MODE == "manual":
        return sorted(MANUAL_ORDERS)

    if ORDER_MODE == "Q3mod4":
        return prime_power_orders(BOUND, odd_only=True, congruence_mod4=3)

    if ORDER_MODE == "basep3mod4":
        orders = []
        for Q in range(3, BOUND):
            data = prime_power_data(Q)
            if data is None:
                continue
            p, _k = data
            if p % 4 == 3:
                orders.append(Q)
        return orders

    raise ValueError(f"Unknown ORDER_MODE={ORDER_MODE!r}")


# ============================================================
# CENTER-ONE LINE HELPERS
# ============================================================


def const(F: FiniteField, c: int) -> int:
    """Return the prime-field constant c inside GF(p^k)."""
    return c % F.p


def lambda_set_for_b(F: FiniteField, b: int):
    one = F.one
    return [
        one,
        F.neg(one),
        b,
        F.neg(b),
        F.add(one, b),
        F.neg(F.add(one, b)),
        F.sub(one, b),
        F.sub(b, one),
    ]


def b_is_valid(F: FiniteField, b: int) -> bool:
    L = lambda_set_for_b(F, b)
    return F.zero not in L and len(set(L)) == 8


def center_one_b3_entries(F: FiniteField, t: int):
    """
    Center-one entries on the line u=t, v=3t.

    The entries are 1 + lambda*t for lambda =
        1, -4, 3, 2, 0, -2, -3, 4, -1.

    Because these lambdas are prime-field constants, scalar multiplication is
    faster than full field multiplication.
    """
    return [
        F.add_one_to_scalar_mul(1, t),
        F.add_one_to_scalar_mul(-4, t),
        F.add_one_to_scalar_mul(3, t),
        F.add_one_to_scalar_mul(2, t),
        F.one,
        F.add_one_to_scalar_mul(-2, t),
        F.add_one_to_scalar_mul(-3, t),
        F.add_one_to_scalar_mul(4, t),
        F.add_one_to_scalar_mul(-1, t),
    ]


def add_line_metadata(w, b=None, t=None):
    if w is not None:
        if b is not None:
            w["b"] = b
        if t is not None:
            w["t"] = t
    return w


def q3_witness_to_text(w):
    F = w["F"]
    lines = []
    if "b" in w:
        lines.append(f"b = {F.elem_str(w['b'])}")
    if "t" in w:
        lines.append(f"t = {F.elem_str(w['t'])}")
    lines.append(witness_to_text(w))
    return "\n".join(lines)


def find_center_one_b3_line_witness(F: FiniteField, n: int, P, root):
    b = const(F, 3)
    if not b_is_valid(F, b):
        return None, "invalid_b"

    for t in range(F.Q):
        entries = center_one_b3_entries(F, t)
        if len(set(entries)) != 9:
            continue
        if all(z in P for z in entries):
            U = t
            V = F.scalar_mul(3, t)
            w = make_witness(F, n, "center_one_line_b3", U, V, entries, root)
            return add_line_metadata(w, b=b, t=t), "found"

    return None, "not_found"


def find_center_one_line_for_b(F: FiniteField, n: int, P, root, b: int):
    if not b_is_valid(F, b):
        return None

    for t in range(F.Q):
        U = t
        V = F.mul(b, t)
        entries = center_one_entries(F, U, V)
        if len(set(entries)) != 9:
            continue
        if all(z in P for z in entries):
            w = make_witness(F, n, "center_one_line_auto_b", U, V, entries, root)
            return add_line_metadata(w, b=b, t=t)

    return None


def candidate_b_values(F: FiniteField):
    """
    Candidate b values for the auto line search.

    Constants are inserted as prime-field constants. For extension fields, alpha
    is also tried early because it often avoids small-characteristic collisions.
    After that, all field elements are tried, so auto_line is exhaustive over
    one-parameter center-one lines u=t, v=b t.
    """
    first = []

    # Constants likely to work in ordinary odd characteristics.
    for c in [3, 4, 5, 6, 2, -2, -3, -4, 7, -7]:
        first.append(const(F, c))

    # For extension fields, try alpha early.
    if F.k > 1:
        first.append(F.alpha)
        first.append(F.add(F.alpha, F.one))
        first.append(F.sub(F.alpha, F.one))

    seen = set()
    for b in first:
        if b not in seen:
            seen.add(b)
            yield b

    for b in range(F.Q):
        if b not in seen:
            seen.add(b)
            yield b


def find_center_one_auto_line_witness(F: FiniteField, n: int, P, root):
    for b in candidate_b_values(F):
        if not b_is_valid(F, b):
            continue
        w = find_center_one_line_for_b(F, n, P, root, b)
        if w is not None:
            return w
    return None


def find_q3_witness(F: FiniteField, mode: str):
    P, root = F.n_power_set_and_roots(EXPONENT_N)

    if mode == "b3_line":
        w, _status = find_center_one_b3_line_witness(F, EXPONENT_N, P, root)
        return w

    if mode == "auto_line":
        return find_center_one_auto_line_witness(F, EXPONENT_N, P, root)

    if mode == "b3_then_auto_line":
        w, _status = find_center_one_b3_line_witness(F, EXPONENT_N, P, root)
        if w is not None:
            return w
        return find_center_one_auto_line_witness(F, EXPONENT_N, P, root)

    if mode == "full_uv":
        return find_center_one_full_witness(F, EXPONENT_N, P, root)

    raise ValueError(f"Unknown MODE={mode!r}")


# ============================================================
# RUNNER
# ============================================================


def run():
    orders = relevant_orders()
    witnesses = []
    failures = []
    errors = []
    start = perf_counter()

    with open(RESULTS_TXT, "w", encoding="utf-8") as out:
        out.write("q ≡ 3 (mod 4) square-field center-one verification / search\n")
        out.write("===========================================================\n")
        out.write(f"bound: Q < {BOUND}\n")
        out.write(f"exponent n: {EXPONENT_N}\n")
        out.write(f"order mode: {ORDER_MODE}\n")
        out.write(f"search mode: {MODE}\n")
        out.write(f"number of fields: {len(orders)}\n")
        out.write(f"orders: {orders}\n\n")
        out.flush()

        for idx, Q in enumerate(orders, start=1):
            t0 = perf_counter()
            out.write("\n" + "="*72 + "\n")
            out.write(f"Checking F_{Q}\n")
            out.flush()

            try:
                F = FiniteField(Q)
                w = find_q3_witness(F, MODE)
                dt = perf_counter() - t0

                if w is None:
                    failures.append(Q)
                    msg = f"F_{Q}: no center-one witness found ({dt:.3f}s)"
                    print(msg)
                    out.write(msg + "\n")
                else:
                    witnesses.append(Q)
                    msg = f"F_{Q}: center-one witness found ({dt:.3f}s), mode={w['mode']}"
                    print(msg)
                    out.write(msg + "\n")
                    out.write(q3_witness_to_text(w) + "\n")

            except Exception as e:
                dt = perf_counter() - t0
                failures.append(Q)
                errors.append((Q, repr(e)))
                msg = f"F_{Q}: ERROR after {dt:.3f}s: {e!r}"
                print(msg)
                out.write(msg + "\n")

            finally:
                out.flush()
                try:
                    del F
                except Exception:
                    pass
                try:
                    del w
                except Exception:
                    pass
                if GC_EVERY and idx % GC_EVERY == 0:
                    gc.collect()

        elapsed = perf_counter() - start
        out.write("\n" + "="*72 + "\n")
        out.write("SUMMARY\n")
        out.write("="*72 + "\n")
        out.write(f"total time: {elapsed:.3f}s\n")
        out.write(f"center-one witnesses found for: {witnesses}\n")
        out.write(f"no center-one witness found: {failures}\n")
        if errors:
            out.write(f"errors: {errors}\n")
        out.write("\nInterpretation note:\n")
        out.write("In line-search modes, a failure only means no line witness was found.\n")
        out.write("For Q ≡ 3 mod 4, rerunning failures with MODE='full_uv' gives a complete\n")
        out.write("center-one check; since -1 is not a square in these fields, center-zero\n")
        out.write("cannot give a distinct square magic square.\n")

    print("\nSUMMARY")
    print("center-one witnesses found for:", witnesses)
    print("no center-one witness found:", failures)
    if errors:
        print("errors:", errors)
    print("results written to:", RESULTS_TXT)


if __name__ == "__main__":
    run()
