from math import isqrt
from time import perf_counter
import gc

# ============================================================
# SETTINGS - edit these and run with:  python center_one_pure_python.py
# ============================================================

BOUND = 553736
ORDER_MODE = "Q3mod4"   # "single", "Q3mod4", or "basep3mod4"
SINGLE_Q = 59
MODE = "b3_then_auto_line"  # "b3_line", "auto_line", "b3_then_auto_line", "full_uv"
PRINT_FIRST_WITNESS = True
# Memory-safe default: store only field orders, not full witness objects.
STORE_FULL_WITNESSES = False
GC_EVERY = 20

# Write progress/results to a text file as the scan runs.
WRITE_RESULTS_TXT = True
RESULTS_TXT = "q3_results.txt"



# ============================================================
# INTEGER / PRIME POWER HELPERS
# ============================================================

def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def factor_int(n):
    fac = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            fac.append((d, e))
        d += 1 if d == 2 else 2
    if n > 1:
        fac.append((n, 1))
    return fac


def prime_power_data(Q):
    fac = factor_int(Q)
    if len(fac) != 1:
        return None
    return fac[0]


def relevant_orders(bound, mode, single_q):
    if mode == "single":
        return [single_q]

    out = []
    for Q in range(2, bound):
        data = prime_power_data(Q)
        if data is None:
            continue
        p, k = data
        if mode == "Q3mod4" and Q % 4 == 3:
            out.append(Q)
        elif mode == "basep3mod4" and p % 4 == 3:
            out.append(Q)
    return out


def prime_divisors(n):
    return [p for p, e in factor_int(n)]


# ============================================================
# POLYNOMIALS OVER F_p, LOW-TO-HIGH COEFFICIENTS
# Used only to find an irreducible polynomial for GF(p^k).
# ============================================================

def poly_trim(f):
    f = list(f)
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    return f


def poly_add(f, g, p):
    n = max(len(f), len(g))
    h = [0] * n
    for i in range(n):
        h[i] = ((f[i] if i < len(f) else 0) + (g[i] if i < len(g) else 0)) % p
    return poly_trim(h)


def poly_sub(f, g, p):
    n = max(len(f), len(g))
    h = [0] * n
    for i in range(n):
        h[i] = ((f[i] if i < len(f) else 0) - (g[i] if i < len(g) else 0)) % p
    return poly_trim(h)


def poly_mul(f, g, p):
    h = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        if a:
            for j, b in enumerate(g):
                if b:
                    h[i+j] = (h[i+j] + a*b) % p
    return poly_trim(h)


def poly_divmod(f, g, p):
    f = poly_trim(f)
    g = poly_trim(g)
    if g == [0]:
        raise ZeroDivisionError
    q = [0] * max(1, (len(f) - len(g) + 1))
    r = f[:]
    inv_lc = pow(g[-1], -1, p)
    while len(r) >= len(g) and r != [0]:
        shift = len(r) - len(g)
        coeff = r[-1] * inv_lc % p
        q[shift] = coeff
        for i in range(len(g)):
            r[shift+i] = (r[shift+i] - coeff*g[i]) % p
        r = poly_trim(r)
    return poly_trim(q), r


def poly_mod(f, m, p):
    return poly_divmod(f, m, p)[1]


def poly_gcd(f, g, p):
    f = poly_trim(f)
    g = poly_trim(g)
    while g != [0]:
        f, g = g, poly_mod(f, g, p)
    if f == [0]:
        return [0]
    inv_lc = pow(f[-1], -1, p)
    return [(c * inv_lc) % p for c in f]


def poly_mul_mod(f, g, m, p):
    return poly_mod(poly_mul(f, g, p), m, p)


def poly_pow_mod(base, exp, m, p):
    result = [1]
    base = poly_mod(base, m, p)
    while exp:
        if exp & 1:
            result = poly_mul_mod(result, base, m, p)
        base = poly_mul_mod(base, base, m, p)
        exp >>= 1
    return result


def is_irreducible_poly(f, p):
    """Rabin irreducibility test for monic f over F_p."""
    f = poly_trim(f)
    n = len(f) - 1
    if n <= 0 or f[-1] % p == 0:
        return False
    if n == 1:
        return True

    x = [0, 1]

    # x^(p^n) == x mod f
    xp_n = poly_pow_mod(x, p**n, f, p)
    if poly_sub(xp_n, x, p) != [0]:
        return False

    # gcd(f, x^(p^(n/r)) - x) = 1 for each prime r | n
    for r in prime_divisors(n):
        xp = poly_pow_mod(x, p**(n // r), f, p)
        h = poly_sub(xp, x, p)
        if poly_gcd(f, h, p) != [1]:
            return False

    return True


def int_to_coeffs(n, p, k):
    coeffs = []
    for _ in range(k):
        coeffs.append(n % p)
        n //= p
    return coeffs


def find_irreducible_poly(p, k):
    """
    Return a monic irreducible polynomial of degree k over F_p,
    represented low-to-high. For k=1, returns x.
    """
    if k == 1:
        return [0, 1]

    # Search monic polynomials c0 + c1*x + ... + c_{k-1}*x^{k-1} + x^k.
    # Skip c0=0 because then x divides the polynomial.
    for n in range(p**k):
        coeffs = int_to_coeffs(n, p, k)
        if coeffs[0] == 0:
            continue
        f = coeffs + [1]
        if is_irreducible_poly(f, p):
            return f

    raise RuntimeError(f"No irreducible polynomial found for p={p}, k={k}")


# ============================================================
# FINITE FIELD GF(p^k), ELEMENTS ENCODED AS INTEGERS 0,...,Q-1
# ============================================================

class FiniteField:
    def __init__(self, Q):
        data = prime_power_data(Q)
        if data is None:
            raise ValueError(f"Q={Q} is not a prime power")
        self.Q = int(Q)
        self.p, self.k = data
        if self.p == 2:
            raise ValueError("This script is for odd finite fields")
        self.powp = [self.p**i for i in range(self.k)]
        # Prime fields are represented by plain integers modulo p.
        # Avoid storing Q singleton coefficient tuples when k = 1.
        if self.k == 1:
            self.modulus = [0, 1]
            self.coeffs = None
        else:
            self.modulus = find_irreducible_poly(self.p, self.k)
            self.coeffs = [tuple(int_to_coeffs(a, self.p, self.k)) for a in range(self.Q)]
        self.zero = 0
        self.one = 1
        self.alpha = self.p if self.k > 1 else 1

    def encode(self, coeffs):
        if self.k == 1:
            return coeffs[0] % self.p if coeffs else 0
        s = 0
        for i, c in enumerate(coeffs[:self.k]):
            s += (c % self.p) * self.powp[i]
        return s

    def add(self, a, b):
        if self.k == 1:
            return (a + b) % self.p
        ca = self.coeffs[a]
        cb = self.coeffs[b]
        return self.encode([(ca[i] + cb[i]) % self.p for i in range(self.k)])

    def neg(self, a):
        if self.k == 1:
            return (-a) % self.p
        ca = self.coeffs[a]
        return self.encode([(-ca[i]) % self.p for i in range(self.k)])

    def sub(self, a, b):
        if self.k == 1:
            return (a - b) % self.p
        ca = self.coeffs[a]
        cb = self.coeffs[b]
        return self.encode([(ca[i] - cb[i]) % self.p for i in range(self.k)])

    def scalar_mul(self, c, a):
        c %= self.p
        if c == 0 or a == 0:
            return 0
        if self.k == 1:
            return (c * a) % self.p
        ca = self.coeffs[a]
        return self.encode([(c * ca[i]) % self.p for i in range(self.k)])

    def add_one_to_scalar_mul(self, c, a):
        c %= self.p
        if self.k == 1:
            return (1 + c * a) % self.p
        ca = self.coeffs[a]
        out = [(c * ca[i]) % self.p for i in range(self.k)]
        out[0] = (out[0] + 1) % self.p
        return self.encode(out)

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        if self.k == 1:
            return (a * b) % self.p

        ca = self.coeffs[a]
        cb = self.coeffs[b]
        prod = [0] * (2*self.k - 1)

        for i, ai in enumerate(ca):
            if ai:
                for j, bj in enumerate(cb):
                    if bj:
                        prod[i+j] = (prod[i+j] + ai*bj) % self.p

        # modulus is m0 + m1*x + ... + m_{k-1}*x^{k-1} + x^k
        m = self.modulus
        for d in range(2*self.k - 2, self.k - 1, -1):
            coeff = prod[d] % self.p
            if coeff:
                shift = d - self.k
                # x^d = x^shift * x^k = -sum m[j] x^(j+shift)
                for j in range(self.k):
                    prod[shift+j] = (prod[shift+j] - coeff*m[j]) % self.p

        return self.encode(prod[:self.k])

    def square_set_and_roots(self):
        S = set()
        root = {}
        if self.k == 1:
            p = self.p
            for a in range(p):
                z = (a * a) % p
                S.add(z)
                if z not in root:
                    root[z] = a
            return S, root

        for a in range(self.Q):
            z = self.mul(a, a)
            S.add(z)
            if z not in root:
                root[z] = a
        return S, root

    def elem_str(self, a):
        if self.k == 1:
            return str(a % self.p)
        cs = self.coeffs[a]
        terms = []
        for i, c in enumerate(cs):
            if c == 0:
                continue
            if i == 0:
                terms.append(str(c))
            elif i == 1:
                if c == 1:
                    terms.append("a")
                else:
                    terms.append(f"{c}*a")
            else:
                if c == 1:
                    terms.append(f"a^{i}")
                else:
                    terms.append(f"{c}*a^{i}")
        return " + ".join(terms) if terms else "0"

    def poly_str(self):
        # modulus low-to-high: m0 + ... + x^k = 0, so alpha satisfies this polynomial.
        terms = []
        for i, c in enumerate(self.modulus):
            c %= self.p
            if c == 0:
                continue
            if i == 0:
                terms.append(str(c))
            elif i == 1:
                terms.append("x" if c == 1 else f"{c}*x")
            else:
                terms.append(f"x^{i}" if c == 1 else f"{c}*x^{i}")
        return " + ".join(terms)


# ============================================================
# CENTER-ONE SEARCHES
# ============================================================

def center_one_entries(F, u, v):
    one = F.one
    return [
        F.add(one, u),          F.sub(F.sub(one, u), v), F.add(one, v),
        F.add(F.sub(one, u), v), one,                    F.sub(F.add(one, u), v),
        F.sub(one, v),          F.add(F.add(one, u), v), F.sub(one, u),
    ]


def print_witness(w):
    F = w["F"]
    print(f"Q = {F.Q} = {F.p}^{F.k}")
    if F.k > 1:
        print(f"Field model: GF({F.Q}) = F_{F.p}[a]/({F.poly_str()})")
    print("mode =", w["mode"])
    if "b" in w:
        print("b =", F.elem_str(w["b"]))
    if "t" in w:
        print("t =", F.elem_str(w["t"]))
    print("u =", F.elem_str(w["u"]))
    print("v =", F.elem_str(w["v"]))

    print("\nMagic square entries:")
    E = w["entries"]
    for i in range(0, 9, 3):
        print([F.elem_str(x) for x in E[i:i+3]])

    print("\nOne choice of square roots:")
    R = w["roots"]
    for i in range(0, 9, 3):
        print([F.elem_str(x) for x in R[i:i+3]])



def witness_to_text(w):
    """Return a text representation of a witness for saving to a file."""
    F = w["F"]
    lines = []
    lines.append(f"F_{F.Q} = GF({F.Q}) = GF({F.p}^{F.k})")
    lines.append(f"characteristic: {F.p}")
    if F.k > 1:
        lines.append(f"field model: GF({F.Q}) = F_{F.p}[a]/({F.poly_str()})")
    lines.append(f"search mode: {w['mode']}")
    if "b" in w:
        lines.append(f"b = {F.elem_str(w['b'])}")
    if "t" in w:
        lines.append(f"t = {F.elem_str(w['t'])}")
    lines.append(f"u = {F.elem_str(w['u'])}")
    lines.append(f"v = {F.elem_str(w['v'])}")

    lines.append("magic square entries:")
    E = w["entries"]
    for i in range(0, 9, 3):
        lines.append("  " + str([F.elem_str(x) for x in E[i:i+3]]))

    lines.append("one choice of square roots:")
    R = w["roots"]
    for i in range(0, 9, 3):
        lines.append("  " + str([F.elem_str(x) for x in R[i:i+3]]))

    return "\n".join(lines)


def write_text_block(handle, title, text=""):
    """Write a visibly separated block and flush immediately."""
    handle.write("\n" + "="*72 + "\n")
    handle.write(title + "\n")
    handle.write("="*72 + "\n")
    if text:
        handle.write(text.rstrip() + "\n")
    handle.flush()


def valid_b(F, b):
    one = F.one
    vals = [
        one,
        F.neg(one),
        b,
        F.neg(b),
        F.add(one, b),
        F.neg(F.add(one, b)),
        F.sub(one, b),
        F.sub(b, one),
    ]
    return 0 not in vals and len(set(vals)) == 8


def find_b3_line_witness(F, S, root):
    # b = 3 in the prime subfield. This is fastest because all multipliers are scalars.
    p = F.p
    scalar_lambdas = [1, -1, 3, -3, 4, -4, -2, 2]
    reduced = [c % p for c in scalar_lambdas]
    if 0 in reduced or len(set(reduced)) != 8:
        return None, "invalid_b3"

    for t in range(1, F.Q):
        entries = [F.add_one_to_scalar_mul(c, t) for c in scalar_lambdas]
        # row-major order corresponding to u=t, v=3t:
        # [1+t, 1-4t, 1+3t, 1+2t, 1, 1-2t, 1-3t, 1+4t, 1-t]
        entries = [entries[0], entries[5], entries[2], entries[7], F.one, entries[6], entries[3], entries[4], entries[1]]

        # For valid b=3 and t != 0, distinctness is automatic, but keep the check.
        if len(set(entries)) != 9:
            continue
        if all(x in S for x in entries):
            return {
                "F": F,
                "mode": "b3_line",
                "b": F.scalar_mul(3, F.one),
                "t": t,
                "u": t,
                "v": F.scalar_mul(3, t),
                "entries": entries,
                "roots": [root[x] for x in entries],
            }, None
    return None, None


def find_auto_line_witness(F, S, root):
    # Try b=3, then alpha, then powers/elements. This is useful in characteristics 3 and 7.
    candidates = []
    b3 = F.scalar_mul(3, F.one)
    candidates.append(b3)
    if F.k > 1:
        candidates.append(F.alpha)
        x = F.alpha
        for _ in range(2, min(F.Q, 30)):
            x = F.mul(x, F.alpha)
            candidates.append(x)
    candidates.extend(range(F.Q))

    seen_b = set()
    one = F.one

    for b in candidates:
        if b in seen_b:
            continue
        seen_b.add(b)
        if not valid_b(F, b):
            continue

        lambdas = [
            one,
            F.neg(one),
            b,
            F.neg(b),
            F.add(one, b),
            F.neg(F.add(one, b)),
            F.sub(one, b),
            F.sub(b, one),
        ]

        for t in range(1, F.Q):
            vals = [F.add(one, F.mul(lam, t)) for lam in lambdas]
            # row-major order for u=t, v=bt:
            # [1+t, 1-(1+b)t, 1+bt, 1+(b-1)t, 1, 1+(1-b)t, 1-bt, 1+(1+b)t, 1-t]
            entries = [vals[0], vals[5], vals[2], vals[7], one, vals[6], vals[3], vals[4], vals[1]]

            if len(set(entries)) != 9:
                continue
            if all(x in S for x in entries):
                return {
                    "F": F,
                    "mode": "auto_line",
                    "b": b,
                    "t": t,
                    "u": t,
                    "v": F.mul(b, t),
                    "entries": entries,
                    "roots": [root[x] for x in entries],
                }
    return None


def find_full_uv_witness(F, S, root):
    one = F.one
    good_u = [u for u in range(F.Q) if F.add(one, u) in S and F.sub(one, u) in S]
    good_v = [v for v in range(F.Q) if F.add(one, v) in S and F.sub(one, v) in S]

    for u in good_u:
        for v in good_v:
            entries = center_one_entries(F, u, v)
            if len(set(entries)) != 9:
                continue
            if all(x in S for x in entries):
                return {
                    "F": F,
                    "mode": "full_uv",
                    "u": u,
                    "v": v,
                    "entries": entries,
                    "roots": [root[x] for x in entries],
                }
    return None


def find_witness_for_Q(Q, mode):
    F = FiniteField(Q)
    S, root = F.square_set_and_roots()

    if mode == "b3_line":
        w, reason = find_b3_line_witness(F, S, root)
        return w, reason

    if mode == "auto_line":
        return find_auto_line_witness(F, S, root), None

    if mode == "b3_then_auto_line":
        w, reason = find_b3_line_witness(F, S, root)
        if w is not None:
            return w, None
        w = find_auto_line_witness(F, S, root)
        return w, reason

    if mode == "full_uv":
        return find_full_uv_witness(F, S, root), None

    raise ValueError("Unknown MODE")


# ============================================================
# RUNNER
# ============================================================

def run():
    orders = relevant_orders(BOUND, ORDER_MODE, SINGLE_Q)
    print("ORDER_MODE =", ORDER_MODE)
    print("MODE =", MODE)
    print("BOUND =", BOUND)
    print("STORE_FULL_WITNESSES =", STORE_FULL_WITNESSES)
    print("WRITE_RESULTS_TXT =", WRITE_RESULTS_TXT)
    if WRITE_RESULTS_TXT:
        print("RESULTS_TXT =", RESULTS_TXT)
    print("orders =", orders[:20], "..." if len(orders) > 20 else "")
    print("number of fields =", len(orders))
    print()

    # Memory-safe scan state.
    # Store only the field orders by default. Full witness dictionaries contain the
    # FiniteField object, square roots, and entries, so retaining many of them can
    # keep large per-field tables alive.
    witness_orders = []
    witness_data = {} if STORE_FULL_WITNESSES else None
    first_witness = None
    failures = []
    reasons = {}
    start_all = perf_counter()

    txt = None
    if WRITE_RESULTS_TXT:
        txt = open(RESULTS_TXT, "w", encoding="utf-8")
        txt.write("Center-one finite-field scan results\n")
        txt.write("====================================\n")
        txt.write(f"ORDER_MODE = {ORDER_MODE}\n")
        txt.write(f"MODE = {MODE}\n")
        txt.write(f"BOUND = {BOUND}\n")
        txt.write(f"SINGLE_Q = {SINGLE_Q}\n")
        txt.write(f"number of fields = {len(orders)}\n")
        txt.write("orders = " + repr(orders) + "\n")
        txt.flush()

    try:
        for idx, Q in enumerate(orders, start=1):
            t0 = perf_counter()
            w = None
            reason = None
            try:
                w, reason = find_witness_for_Q(Q, MODE)
            except Exception as e:
                dt = perf_counter() - t0
                msg = f"F_{Q}: ERROR {e} ({dt:.3f}s)"
                print(msg)
                failures.append(Q)
                reasons[Q] = repr(e)
                if txt is not None:
                    write_text_block(txt, f"F_{Q}: ERROR", f"time: {dt:.3f}s\nerror: {repr(e)}")
                if GC_EVERY and idx % GC_EVERY == 0:
                    gc.collect()
                continue

            dt = perf_counter() - t0
            if w is None:
                msg = f"F_{Q}: no center-one witness found ({dt:.3f}s)" + (f" [{reason}]" if reason else "")
                print(msg)
                failures.append(Q)
                if reason:
                    reasons[Q] = reason
                if txt is not None:
                    body = f"time: {dt:.3f}s"
                    if reason:
                        body += f"\nreason: {reason}"
                    write_text_block(txt, f"F_{Q}: NO CENTER-ONE WITNESS FOUND", body)
            else:
                msg = f"F_{Q}: witness found ({dt:.3f}s), mode={w['mode']}"
                print(msg)
                witness_orders.append(Q)

                # Write the full witness immediately before releasing memory.
                if txt is not None:
                    body = f"time: {dt:.3f}s\n" + witness_to_text(w)
                    write_text_block(txt, f"F_{Q}: WITNESS FOUND", body)

                if STORE_FULL_WITNESSES:
                    witness_data[Q] = w
                elif PRINT_FIRST_WITNESS and first_witness is None:
                    # Keep only one full witness so print_witness still works.
                    first_witness = w
                else:
                    # Drop the witness immediately. This releases its FiniteField object.
                    del w

            # Encourage prompt release of large sets/dicts created inside find_witness_for_Q.
            if GC_EVERY and idx % GC_EVERY == 0:
                gc.collect()

    finally:
        if txt is not None:
            elapsed = perf_counter() - start_all
            write_text_block(
                txt,
                "SUMMARY",
                "total time: %.3fs\n" % elapsed
                + "witnesses found:\n"
                + repr(witness_orders)
                + "\n\nno center-one witness found:\n"
                + repr(failures)
                + ("\n\nnotes:\n" + "\n".join(f"F_{Q}: {reasons[Q]}" for Q in sorted(reasons)) if reasons else "")
            )
            txt.close()

    print("\n========================================")
    print("SUMMARY")
    print("========================================")
    print("total time: %.3fs" % (perf_counter() - start_all))
    print("witnesses found:")
    print(witness_orders)
    print("\nno center-one witness found:")
    print(failures)
    if reasons:
        print("\nnotes:")
        for Q in sorted(reasons):
            print(f"F_{Q}: {reasons[Q]}")

    if WRITE_RESULTS_TXT:
        print(f"\nResults written to: {RESULTS_TXT}")

    if PRINT_FIRST_WITNESS:
        if STORE_FULL_WITNESSES and witness_data:
            first_Q = witness_orders[0]
            first_witness_to_print = witness_data[first_Q]
        else:
            first_witness_to_print = first_witness

        if first_witness_to_print is not None:
            print("\n========================================")
            print("FIRST WITNESS")
            print("========================================")
            print_witness(first_witness_to_print)


if __name__ == "__main__":
    run()
