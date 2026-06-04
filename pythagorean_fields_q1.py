"""
Small q ≡ 1 (mod 4) square-field verification.

This script checks every odd prime power Q < 77 with Q ≡ 1 (mod 4).
The paper's character-sum bound proves all larger Q ≡ 1 (mod 4) fields are
non-Parker. The normalized search below is complete: center 0 plus center 1
covers all magic squares of square entries up to scaling by the center.

Run:
    python pythagorean_fields_q1.py

Output:
    q1_results.txt
"""
from time import perf_counter
import gc

from pythagorean_fields_common import (
    FiniteField,
    prime_power_orders,
    find_normalized_witness,
    witness_to_text,
)

BOUND = 77
EXPONENT_N = 2
RESULTS_TXT = "q1_results.txt"
GC_EVERY = 5


def run():
    orders = prime_power_orders(BOUND, odd_only=True, congruence_mod4=1)
    witnesses = []
    parker = []
    start = perf_counter()

    with open(RESULTS_TXT, "w", encoding="utf-8") as out:
        out.write("q ≡ 1 (mod 4) square-field verification\n")
        out.write("==========================================\n")
        out.write(f"bound: Q < {BOUND}\n")
        out.write(f"exponent n: {EXPONENT_N}\n")
        out.write(f"orders: {orders}\n\n")
        out.flush()

        for idx, Q in enumerate(orders, start=1):
            t0 = perf_counter()
            F = FiniteField(Q)
            w = find_normalized_witness(F, EXPONENT_N)
            dt = perf_counter() - t0

            out.write("\n" + "="*72 + "\n")
            if w is None:
                parker.append(Q)
                msg = f"F_{Q}: PARKER / no normalized square witness ({dt:.3f}s)"
                print(msg)
                out.write(msg + "\n")
            else:
                witnesses.append(Q)
                msg = f"F_{Q}: non-Parker witness found ({dt:.3f}s), mode={w['mode']}"
                print(msg)
                out.write(msg + "\n")
                out.write(witness_to_text(w) + "\n")
            out.flush()

            del F, w
            if GC_EVERY and idx % GC_EVERY == 0:
                gc.collect()

        elapsed = perf_counter() - start
        out.write("\n" + "="*72 + "\n")
        out.write("SUMMARY\n")
        out.write("="*72 + "\n")
        out.write(f"total time: {elapsed:.3f}s\n")
        out.write(f"non-Parker fields: {witnesses}\n")
        out.write(f"Parker fields: {parker}\n")

    print("\nSUMMARY")
    print("non-Parker fields:", witnesses)
    print("Parker fields:", parker)
    print("results written to:", RESULTS_TXT)


if __name__ == "__main__":
    run()
