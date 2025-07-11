import itertools

def get_exc(perm):
    """Calculates the number of excedances in a permutation."""
    # A permutation (a_1, ..., a_n) is given as a tuple.
    # An excedance at index i is where a_{i+1} > i+1.
    return sum(1 for i, p_i in enumerate(perm) if p_i > i + 1)

def is_derangement(perm):
    """Checks if a permutation is a derangement."""
    # A derangement has no fixed points, i.e., a_{i+1} != i+1.
    return all(p_i != i + 1 for i, p_i in enumerate(perm))

def poly_to_str(coeffs, var='t'):
    """Converts a dictionary of coefficients to a string."""
    if not coeffs:
        return "0"
    terms = []
    for power in sorted(coeffs.keys()):
        coeff = coeffs[power]
        if power == 0:
            terms.append(str(coeff))
        elif power == 1:
            terms.append(f"{coeff}{var}" if coeff > 1 else var)
        else:
            terms.append(f"{coeff}{var}^{power}" if coeff > 1 else f"{var}^{power}")
    return " + ".join(reversed(terms)) # conventional descending powers

print("Step-by-step analysis:\n")

# --- Part (a) ---
print("(a) Analysis of the identity H(U_{n-1, E})(t) = t^(n-1) * d_n(t)")
n_test = 3
print(f"We test the identity for n = {n_test}.")

# Calculate H(U_{n-1, E})(t) = A_{n-1}(t)
size_A = n_test - 1
A_coeffs = {}
for p in itertools.permutations(range(1, size_A + 1)):
    exc = get_exc(p)
    A_coeffs[exc] = A_coeffs.get(exc, 0) + 1
A_poly_str = poly_to_str(A_coeffs)
print(f"The left side H(U_{{n_test-1}}, E)(t) = A_{n_test-1}(t) = {A_poly_str}")

# Calculate d_n(t)
size_d = n_test
d_coeffs = {}
for p in itertools.permutations(range(1, size_d + 1)):
    if is_derangement(p):
        exc = get_exc(p)
        d_coeffs[exc] = d_coeffs.get(exc, 0) + 1
d_poly_str = poly_to_str(d_coeffs)
print(f"The derangement polynomial d_{n_test}(t) = {d_poly_str}")

# Calculate the right side of the identity: t^(n-1) * d_n(t)
rhs_poly_str = f"t^{n_test-1} * ({d_poly_str})"
print(f"The right side is {rhs_poly_str} = t^2 * (t^2 + t) = t^4 + t^3")

answer_a_identity = "No"
print(f"\nSince {A_poly_str} != t^4 + t^3, the identity is false.")
answer_a_degree = "n-2"
print(f"The degree of H(U_{n-1, E})(t) = A_{n-1}(t) is the maximum number of excedances in a permutation in S_{n-1}, which is n-2.\n")

# --- Part (b) ---
print("(b) Analysis of the leading coefficient of d_n(t)")
print("The leading term of d_n(t) corresponds to the derangement(s) with the maximum number of excedances.")
print("The maximum number of excedances in S_n is n-1. The unique permutation with n-1 excedances is sigma = (2, 3, ..., n, 1).")
print("This permutation is a derangement for all n >= 2, as sigma(i) = i+1 for i < n, and sigma(n) = 1.")
print("Since there is exactly one such derangement, the leading coefficient is 1.")
answer_b = "Yes"
print("The statement is true.\n")


# --- Part (c) ---
print("(c) Calculation of d_3(1)")
# d_coeffs for n=3 was already computed.
# d_3(t) = t + t^2
val_c = 0
calc_str_terms = []
for power, coeff in sorted(d_coeffs.items()):
    val_c += coeff
    calc_str_terms.append(f"{coeff}")
calc_str = " + ".join(calc_str_terms)

# This fulfills the "output each number in the final equation" requirement
final_calc_str = f"d_3(1) is the sum of the coefficients of d_3(t) = {d_poly_str}.\nSo, d_3(1) = {calc_str} = {val_c}"
print(final_calc_str)
answer_c = val_c

# --- Final Answer ---
final_answer_string = f"\nFinal Answer: (a) {answer_a_identity}, {answer_a_degree}; (b) {answer_b}; (c) {answer_c}"
print(final_answer_string)
print("\n<<<" + f"(a) {answer_a_identity}, {answer_a_degree}; (b) {answer_b}; (c) {answer_c}" + ">>>")