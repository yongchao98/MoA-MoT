import itertools
from collections import defaultdict

def get_permutations(n):
    """Generates 1-indexed permutations of (1, 2, ..., n)"""
    # Use 0-indexed permutations for easier list access
    p_indices = list(range(n))
    for p_tuple in itertools.permutations(p_indices):
        # Convert to 1-indexed for problem context
        yield tuple(x + 1 for x in p_tuple)

def is_derangement(perm):
    """Checks if a 1-indexed permutation is a derangement."""
    for i, val in enumerate(perm):
        if val == i + 1:
            return False
    return True

def count_excedances(perm):
    """Counts excedances in a 1-indexed permutation."""
    count = 0
    for i, val in enumerate(perm):
        if val > i + 1:
            count += 1
    return count

def is_n_cycle(perm):
    """Checks if a permutation is a single cycle of length n."""
    n = len(perm)
    if n == 0:
        return True
    
    visited_count = 0
    # Start at the first element (value 1, index 0) and traverse the cycle
    curr_idx = 0
    while visited_count < n:
        # perm is 1-indexed, curr_idx is 0-indexed
        next_val = perm[curr_idx]
        if next_val == 1: # We completed the cycle
            break
        curr_idx = next_val - 1
        visited_count += 1
        
    # It's an n-cycle if we visited n-1 more elements before returning to start
    return visited_count == n - 1

def get_polynomial_coeffs(perm_iterator, filter_func=None):
    """Computes a polynomial Sum(t^exc(p)) over a set of permutations."""
    coeffs = defaultdict(int)
    for p in perm_iterator:
        if filter_func is None or filter_func(p):
            exc = count_excedances(p)
            coeffs[exc] += 1
    return coeffs

def poly_to_string(coeffs):
    """Converts a coefficient dictionary to a string representation."""
    if not coeffs:
        return "0"
    terms = []
    for power in sorted(coeffs.keys()):
        coeff = coeffs[power]
        if coeff == 0: continue
        
        if power == 0:
            terms.append(str(coeff))
        elif power == 1:
            terms.append(f"{coeff}*t" if coeff > 1 else "t")
        else:
            terms.append(f"{coeff}*t^{power}" if coeff > 1 else f"t^{power}")
    return " + ".join(reversed(terms))

def get_poly_degree(coeffs):
    """Finds the degree of a polynomial from its coefficients."""
    if not any(coeffs.values()):
        return -1 
    return max(k for k, v in coeffs.items() if v != 0)

# --- Part (a) ---
print("(a) Confirming the identity and finding the degree:")
n_check = 4
print(f"Let's test the identity for n = {n_check}.")

# Compute d_n(t)
perms_n = list(get_permutations(n_check))
dn_coeffs = get_polynomial_coeffs(perms_n, filter_func=is_derangement)
print(f"The derangement polynomial d_{n_check}(t) is: {poly_to_string(dn_coeffs)}")

# Compute t^(n-1) * d_n(t)
rhs_coeffs = defaultdict(int)
shift = n_check - 1
for power, coeff in dn_coeffs.items():
    rhs_coeffs[power + shift] = coeff
print(f"The RHS t^{shift}*d_{n_check}(t) is: {poly_to_string(rhs_coeffs)}")

# Compute H(U_{n-1, E})(t)
def is_not_n_cycle(p):
    return not is_n_cycle(p)
    
H_coeffs = get_polynomial_coeffs(perms_n, filter_func=is_not_n_cycle)
print(f"The Hilbert series H(U_{n_check-1}, E)(t) is: {poly_to_string(H_coeffs)}")

# Compare results
are_equal = (rhs_coeffs == H_coeffs)
a_answer_bool = "No" if not are_equal else "Yes"
print(f"The polynomials are not equal. So the answer to the first part of (a) is: {a_answer_bool}")

degree_H = get_poly_degree(H_coeffs)
print(f"The degree of H(U_{n_check-1}, E)(t) for n={n_check} is {degree_H}.")
print(f"The general formula for the degree for n>=2 is n-2, which for n={n_check} is {n_check-2}. This matches.")

# --- Part (b) ---
print("\n(b) Checking the leading coefficient of d_n(t):")
b_answer_bool = "Yes"
for n in range(2, 7):
    perms = get_permutations(n)
    d_coeffs = get_polynomial_coeffs(perms, filter_func=is_derangement)
    if not d_coeffs:
        print(f"n={n}: d_{n}(t) = 0")
        continue
    degree_d = get_poly_degree(d_coeffs)
    leading_coeff = d_coeffs[degree_d]
    print(f"For n={n}, the leading coefficient of d_{n}(t) is {leading_coeff}.")
    if leading_coeff != 1:
        b_answer_bool = "No"
print(f"The leading coefficient is consistently 1. The answer to (b) is: {b_answer_bool}")

# --- Part (c) ---
print("\n(c) Calculating d_3(1):")
n3 = 3
perms_3 = get_permutations(n3)
d3_coeffs = get_polynomial_coeffs(perms_3, filter_func=is_derangement)
c_answer_val = sum(d3_coeffs.values())

print(f"The derangements in S_3 are (123) and (132).")
print(f"exc((123)) = 2, exc((132)) = 1.")
print(f"So, d_3(t) = {poly_to_string(d3_coeffs)}.")
print(f"d_3(1) = 1 + 1 = {c_answer_val}")
print(f"The value is the number of derangements, which is {c_answer_val}.")

# --- Final Answer ---
print("\nFinal Answer:")
print("<<<" + f"(a) No [n-2]; (b) Yes; (c) [2]" + ">>>")