"""
Dependency-free finite-field utilities for magic-square searches.

Finite fields GF(p^k) are represented by integers 0,...,p^k-1, interpreted as
polynomials in alpha with coefficients in F_p. For k>1, alpha is a root of the
first monic irreducible polynomial found over F_p.
"""
from math import gcd


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


def prime_divisors(n):
    return [p for p, _ in factor_int(n)]


def prime_power_orders(bound, *, odd_only=True, congruence_mod4=None):
    out = []
    for Q in range(2, bound):
        data = prime_power_data(Q)
        if data is None:
            continue
        if odd_only and Q % 2 == 0:
            continue
        if congruence_mod4 is not None and Q % 4 != congruence_mod4:
            continue
        out.append(Q)
    return out


# ============================================================
# Polynomials over F_p, coefficients low-to-high
# ============================================================


def poly_trim(f):
    f = list(f)
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    return f


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
    q = [0] * max(1, len(f) - len(g) + 1)
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
    if poly_sub(poly_pow_mod(x, p**n, f, p), x, p) != [0]:
        return False
    for r in prime_divisors(n):
        h = poly_sub(poly_pow_mod(x, p**(n // r), f, p), x, p)
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
    if k == 1:
        return [0, 1]
    for n in range(p**k):
        coeffs = int_to_coeffs(n, p, k)
        if coeffs[0] == 0:
            continue
        f = coeffs + [1]
        if is_irreducible_poly(f, p):
            return f
    raise RuntimeError(f"No irreducible polynomial found for p={p}, k={k}")


# ============================================================
# Finite field GF(p^k)
# ============================================================


class FiniteField:
    def __init__(self, Q):
        data = prime_power_data(Q)
        if data is None:
            raise ValueError(f"Q={Q} is not a prime power")
        self.Q = int(Q)
        self.p, self.k = data
        if self.p == 2:
            raise ValueError("This helper is intended for odd finite fields")
        self.powp = [self.p**i for i in range(self.k)]
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
            return (1 + c*a) % self.p
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
        m = self.modulus
        for d in range(2*self.k - 2, self.k - 1, -1):
            coeff = prod[d] % self.p
            if coeff:
                shift = d - self.k
                for j in range(self.k):
                    prod[shift+j] = (prod[shift+j] - coeff*m[j]) % self.p
        return self.encode(prod[:self.k])

    def pow(self, a, e):
        result = 1
        base = a
        while e:
            if e & 1:
                result = self.mul(result, base)
            base = self.mul(base, base)
            e >>= 1
        return result

    def n_power_set_and_roots(self, n):
        P = set()
        root = {}
        for a in range(self.Q):
            z = self.pow(a, n)
            P.add(z)
            if z not in root:
                root[z] = a
        return P, root

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
                terms.append("a" if c == 1 else f"{c}*a")
            else:
                terms.append(f"a^{i}" if c == 1 else f"{c}*a^{i}")
        return " + ".join(terms) if terms else "0"

    def poly_str(self):
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
# Magic-square entry formulas
# ============================================================


def center_zero_entries(F, U, V):
    return [
        U,                         F.neg(F.add(U, V)), V,
        F.sub(V, U),               0,                  F.sub(U, V),
        F.neg(V),                  F.add(U, V),         F.neg(U),
    ]


def center_one_entries(F, U, V):
    one = F.one
    return [
        F.add(one, U),             F.sub(F.sub(one, U), V), F.add(one, V),
        F.add(F.sub(one, U), V),   one,                      F.sub(F.add(one, U), V),
        F.sub(one, V),             F.add(F.add(one, U), V),   F.sub(one, U),
    ]


def witness_to_text(w):
    F = w["F"]
    lines = []
    lines.append(f"F_{F.Q} = GF({F.Q}) = GF({F.p}^{F.k})")
    lines.append(f"characteristic: {F.p}")
    if F.k > 1:
        lines.append(f"field model: GF({F.Q}) = F_{F.p}[a]/({F.poly_str()})")
    lines.append(f"exponent n: {w['n']}")
    lines.append(f"search mode: {w['mode']}")
    if "x" in w:
        lines.append(f"x = {F.elem_str(w['x'])}")
    if "U" in w:
        lines.append(f"U = {F.elem_str(w['U'])}")
    if "V" in w:
        lines.append(f"V = {F.elem_str(w['V'])}")
    lines.append("magic square entries:")
    E = w["entries"]
    for i in range(0, 9, 3):
        lines.append("  " + str([F.elem_str(x) for x in E[i:i+3]]))
    lines.append(f"one choice of {w['n']}-th roots:")
    R = w["roots"]
    for i in range(0, 9, 3):
        lines.append("  " + str([F.elem_str(x) for x in R[i:i+3]]))
    return "\n".join(lines)


def make_witness(F, n, mode, U, V, entries, root, x=None):
    w = {
        "F": F,
        "n": n,
        "mode": mode,
        "U": U,
        "V": V,
        "entries": entries,
        "roots": [root[z] for z in entries],
    }
    if x is not None:
        w["x"] = x
    return w


def find_center_zero_line_witness(F, n, P, root):
    # M(0,x,1), the one-parameter center-zero family used in the paper.
    one = F.one
    for x in range(F.Q):
        entries = center_zero_entries(F, x, one)
        if len(set(entries)) == 9 and all(z in P for z in entries):
            return make_witness(F, n, "center_zero_line_M0(x)", x, one, entries, root, x=x)
    return None


def find_center_zero_full_witness(F, n, P, root):
    good = [a for a in range(F.Q) if a in P]
    for U in good:
        for V in good:
            entries = center_zero_entries(F, U, V)
            if len(set(entries)) == 9 and all(z in P for z in entries):
                return make_witness(F, n, "center_zero_full_UV", U, V, entries, root)
    return None


def find_center_one_full_witness(F, n, P, root):
    one = F.one
    good_u = [u for u in range(F.Q) if F.add(one, u) in P and F.sub(one, u) in P]
    good_v = [v for v in range(F.Q) if F.add(one, v) in P and F.sub(one, v) in P]
    for U in good_u:
        for V in good_v:
            entries = center_one_entries(F, U, V)
            if len(set(entries)) == 9 and all(z in P for z in entries):
                return make_witness(F, n, "center_one_full_UV", U, V, entries, root)
    return None


def find_normalized_witness(F, n):
    """
    Complete normalized search for n-th powers.

    If a magic square of n-th powers has center E=0, it appears in center-zero form.
    If E is nonzero, then E is itself an n-th power, so scaling by E^{-1} preserves
    n-th powers and normalizes the center to 1. Thus center-zero plus center-one is
    complete for existence/nonexistence.
    """
    P, root = F.n_power_set_and_roots(n)

    # Fast first pass. In the square q≡1 mod 4 case and in the cube case this often
    # finds the expected witness immediately.
    w = find_center_zero_line_witness(F, n, P, root)
    if w is not None:
        return w

    w = find_center_zero_full_witness(F, n, P, root)
    if w is not None:
        return w

    return find_center_one_full_witness(F, n, P, root)
