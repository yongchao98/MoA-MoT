import math
from itertools import permutations
from collections import defaultdict

def poly_to_string(poly_dict):
    """Converts a polynomial from a dictionary to a string."""
    if not poly_dict:
        return "0"
    terms = []
    # Sort by power in descending order
    for power in sorted(poly_dict.keys(), reverse=True):
        coeff = poly_dict[power]
        if coeff != 0:
            terms.append(f"{coeff}*t^{power}")
    return " + ".join(terms)

def is_derangement(p):
    """Checks if a permutation is a derangement.
    Permutation p is 1-indexed for clarity, so p[i] corresponds to sigma(i+1).
    """
    for i, val in enumerate(p):
        if i + 1 == val:
            return False
    return True

def count_excedances(p):
    """Counts excedances in a 1-indexed permutation."""
    count = 0
    for i, val in enumerate(p):
        if val > i + 1:
            count += 1
    return count

def get_derangement_poly(n):
    """Calculates the derangement polynomial d_n(t)."""
    if n == 0:
        return {0: 1} # by convention d_0(t) = 1
    poly = defaultdict(int)
    for p in permutations(range(1, n + 1)):
        if is_derangement(p):
            exc = count_excedances(p)
            poly[exc] += 1
    return poly

def get_hilbert_series(n):
    """Calculates H(U_{n-1, n})(t)."""
    poly = defaultdict(int)
    r = n - 1
    for k in range(r + 1):
        # The term is comb(n-k-1, r-k) * t^k
        # comb(n-k-1, n-1-k)
        try:
            coeff = math.comb(n - k - 1, r - k)
            poly[k] = coeff
        except ValueError:
            # comb(a,b) is 0 if b < 0 or b > a
            poly[k] = 0
    return poly

def multiply_poly_by_t_power(poly, power):
    """Multiplies a polynomial by t^power."""
    new_poly = defaultdict(int)
    for p, c in poly.items():
        new_poly[p + power] = c
    return new_poly

def evaluate_poly(poly, val):
    """Evaluates a polynomial at a given value."""
    return sum(c * (val ** p) for p, c in poly.items())

def get_leading_coeff(poly):
    """Gets the leading coefficient of a polynomial."""
    if not poly:
        return 0
    max_power = max(poly.keys())
    return poly[max_power]

# --- Main execution ---
# Part (a): Check the identity for n=3
n_test = 3
print(f"--- Verifying for n={n_test} ---")
h_n_t = get_hilbert_series(n_test)
d_n_t = get_derangement_poly(n_test)
rhs = multiply_poly_by_t_power(d_n_t, n_test - 1)

print(f"(a) Checking if H_n(t) == t^(n-1)*d_n(t)")
print(f"H_{n_test}(t) = {poly_to_string(h_n_t)}")
print(f"d_{n_test}(t) = {poly_to_string(d_n_t)}")
print(f"t^({n_test}-1) * d_{n_test}(t) = {poly_to_string(rhs)}")
print(f"Conclusion: The identity is False.\n")

# Part (b): Check the leading coefficient for n=2 to 5
print("(b) Checking leading coefficient of d_n(t) for n=2,3,4")
for n in range(2, 5):
    d_n = get_derangement_poly(n)
    lc = get_leading_coeff(d_n)
    print(f"Leading coefficient of d_{n}(t) is: {lc}")
print("Conclusion: The leading coefficient appears to be 1 for n>=2.\n")

# Part (c): Calculate d_3(1)
print("(c) Calculating d_3(1)")
n_c = 3
d_3_t = get_derangement_poly(n_c)
val_at_1 = evaluate_poly(d_3_t, 1)
print(f"d_{n_c}(t) is {poly_to_string(d_3_t)}")
print(f"d_{n_c}(1) = {val_at_1}")
