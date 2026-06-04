"""
Odd-characteristic 3-Parker/cube-field verification below the theoretical bound.

For n=3, the paper proves every odd finite field F_Q with Q >= 1037 contains a
3x3 magic square of distinct cubes. Therefore it remains to check odd prime
powers Q < 1037. This script performs a complete normalized search for those Q.

The normalized search is complete: if the center E is 0, the square is in the
center-zero family; if E is nonzero, then E is itself a cube, so scaling by E^{-1}
preserves cubes and normalizes the center to 1.

Run:
    python pythagorean_fields_n3.py

Output:
    n3_results.txt
"""
from time import perf_counter
import gc

from pythagorean_fields_common import (
    FiniteField,
    prime_power_orders,
    find_normalized_witness,
    witness_to_text,
)

BOUND = 1037
EXPONENT_N = 3
RESULTS_TXT = "n3_results.txt"
GC_EVERY = 10


def run():
    orders = prime_power_orders(BOUND, odd_only=True, congruence_mod4=None)
    witnesses = []
    parker = []
    start = perf_counter()

    with open(RESULTS_TXT, "w", encoding="utf-8") as out:
        out.write("Odd-characteristic 3-Parker / cube-field verification\n")
        out.write("======================================================\n")
        out.write(f"bound: Q < {BOUND}\n")
        out.write(f"exponent n: {EXPONENT_N}\n")
        out.write(f"number of fields: {len(orders)}\n")
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
                msg = f"F_{Q}: 3-PARKER / no normalized cube witness ({dt:.3f}s)"
                print(msg)
                out.write(msg + "\n")
            else:
                witnesses.append(Q)
                msg = f"F_{Q}: non-3-Parker witness found ({dt:.3f}s), mode={w['mode']}"
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
        out.write(f"non-3-Parker odd fields: {witnesses}\n")
        out.write(f"3-Parker odd fields: {parker}\n")
        out.write("\nNote: characteristic two is excluded here; in characteristic two every 3x3\n")
        out.write("magic square has repeated entries, so the distinct-entry problem has permanent\n")
        out.write("obstructions there.\n")

    print("\nSUMMARY")
    print("non-3-Parker odd fields:", witnesses)
    print("3-Parker odd fields:", parker)
    print("results written to:", RESULTS_TXT)


if __name__ == "__main__":
    run()
