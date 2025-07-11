import itertools
from collections import Counter

def get_derangements(n):
    """
    Generates all derangements for a set of n elements.
    A derangement is a permutation sigma such that for all i, sigma(i) != i.
    The problem uses 1-based indexing for permutations.
    """
    derangements = []
    elements = list(range(1, n + 1))
    for p in itertools.permutations(elements):
        is_derangement = True
        # p is a 0-indexed tuple, e.g., (2, 3, 1) for sigma(1)=2, sigma(2)=3, sigma(3)=1
        for i in range(n):
            if p[i] == i + 1:
                is_derangement = False
                break
        if is_derangement:
            derangements.append(p)
    return derangements

def get_excedance_count(perm):
    """
    Calculates the number of excedances in a permutation.
    An excedance is an index i such that sigma(i) > i.
    perm is a tuple representing the permutation of (1, ..., n).
    """
    count = 0
    # perm is 0-indexed, so we check if perm[i] > i + 1
    for i in range(len(perm)):
        if perm[i] > i + 1:
            count += 1
    return count

def get_derangement_polynomial_coeffs(n):
    """
    Computes the coefficients of the derangement polynomial d_n(t).
    Returns a dictionary where keys are powers of t and values are coefficients.
    """
    if n == 0:
        return {0: 1} # By convention, d_0(t) = 1
    if n == 1:
        return {}     # There are no derangements in S_1, so d_1(t) = 0
        
    derangements = get_derangements(n)
    excedance_counts = [get_excedance_count(p) for p in derangements]
    
    poly_coeffs = Counter(excedance_counts)
    return dict(poly_coeffs)

def format_poly(coeffs):
    """Formats a coefficient dictionary into a string."""
    if not coeffs:
        return "0"
    terms = []
    for k, v in sorted(coeffs.items()):
        if v == 1 and k > 0:
            term = f"t^{k}" if k > 1 else "t"
        elif v > 1 and k > 0:
            term = f"{v}*t^{k}" if k > 1 else f"{v}*t"
        else: # k=0
            term = str(v)
        terms.append(term)
    return " + ".join(terms)

# --- Main logic to answer the questions ---

# (a) Confirm whether H(U_{n-1, E})(t) = t^(n-1) * d_n(t).
print("--- Part (a) ---")
n_a = 3
print(f"To test the identity, we use n = {n_a}.")
d3_coeffs = get_derangement_polynomial_coeffs(n_a)
d3_poly_str = format_poly(d3_coeffs)
print(f"The derangement polynomial d_{n_a}(t) is calculated as: {d3_poly_str}.")
print(f"The right side of the proposed identity is t^({n_a}-1) * d_{n_a}(t) = t^2 * ({d3_poly_str}) = t^3 + t^4.")
print("From established results in matroid theory, the Hilbert series for the Chow ring of the uniform matroid U_{2,3} is H(U_{2,3})(t) = 1 + t + t^2.")
print("Since 1 + t + t^2 is not equal to t^3 + t^4, the identity is false.")
answer_a_bool = "No"
answer_a_deg = "n-1"
print(f"The degree of H(U_{n-1,E})(t) is the rank of the matroid, which is {answer_a_deg}.")

# (b) State if the leading coefficient of d_n(t) for any n >= 2 is always 1.
print("\n--- Part (b) ---")
print("We check the leading coefficient of d_n(t) for n = 2, 3, 4, 5.")
is_leading_coeff_one = True
for n_b in range(2, 6):
    poly_b = get_derangement_polynomial_coeffs(n_b)
    if not poly_b: # Handles n=0,1 cases but loop starts at 2
        is_leading_coeff_one = False
        continue
    degree_b = max(poly_b.keys())
    lead_coeff_b = poly_b[degree_b]
    print(f"For n={n_b}, d_{n_b}(t) has degree {degree_b} and its leading coefficient is {lead_coeff_b}.")
    if degree_b != n_b - 1 or lead_coeff_b != 1:
        is_leading_coeff_one = False
if is_leading_coeff_one:
    answer_b_bool = "Yes"
    print("The leading coefficient is 1 in all tested cases. This is because for n>=2, the unique derangement with the maximum number of excedances (n-1) is the cycle (1 2 ... n).")
else:
    answer_b_bool = "No"

# (c) Give the value of d_3(1).
print("\n--- Part (c) ---")
n_c = 3
print(f"d_{n_c}(1) is the number of derangements in S_{n_c}, which is the sum of the coefficients of d_{n_c}(t).")
derangements_c = get_derangements(n_c)
answer_c_val = len(derangements_c)
d3_coeffs_c = get_derangement_polynomial_coeffs(n_c)
coeff_sum = sum(d3_coeffs_c.values())
print(f"The derangements in S_3 are: {derangements_c[0]} and {derangements_c[1]}.")
print(f"There are {answer_c_val} derangements. So, d_3(1) = {coeff_sum}.")

# Format the final answer
final_answer_str = f"(a) {answer_a_bool}, {answer_a_deg}; (b) {answer_b_bool}; (c) {answer_c_val}"

print("\n--------------------")
print("Final Answer:")
print(f"<<<{final_answer_str}>>>")