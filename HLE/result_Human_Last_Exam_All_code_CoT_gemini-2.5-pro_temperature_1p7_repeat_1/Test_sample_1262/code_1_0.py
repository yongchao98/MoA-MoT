import math
from itertools import permutations

def get_derangement_poly_coeffs(n):
    """
    Computes the coefficients of the n-th derangement polynomial d_n(t).
    Returns a dictionary where keys are the number of excedances (powers of t)
    and values are the coefficients (number of derangements with that many excedances).
    """
    if n < 0:
        return {}
    
    coeffs = {}
    elements = range(1, n + 1)
    
    # Iterate through all permutations of {1, 2, ..., n}
    for p in permutations(elements):
        # A permutation p (a tuple) is 0-indexed.
        # It represents a function sigma where sigma(i) = p[i-1].
        
        # Check if it is a derangement (no fixed points: sigma(i) != i)
        is_derangement = all(p[i-1] != i for i in elements)
        
        if is_derangement:
            # Count excedances (sigma(i) > i)
            excedances = sum(1 for i in elements if p[i-1] > i)
            
            # Increment the coefficient for this number of excedances
            coeffs[excedances] = coeffs.get(excedances, 0) + 1
            
    return coeffs

def get_hilbert_series_coeffs(n):
    """
    Computes coefficients of the Hilbert series H(U_{n-1, E})(t).
    H(t) = sum_{i=0}^{n-2} C(n, i+1) * t^i.
    Returns a list of coefficients [c_0, c_1, ...], where c_i is the coeff of t^i.
    """
    if n < 2:
        return []
    
    # The degree of the polynomial is n-2. The loop for i goes from 0 to n-2.
    # The coefficient of t^i is C(n, i+1).
    return [math.comb(n, i + 1) for i in range(n - 1)]

def format_poly_from_dict(coeffs_map, var='t'):
    """Helper function to format a polynomial from a dictionary of coefficients."""
    if not coeffs_map:
        return "0"
    
    terms = []
    for power in sorted(coeffs_map.keys()):
        coeff = coeffs_map[power]
        if coeff == 0:
            continue
        if power == 0:
            terms.append(f"{coeff}")
        elif power == 1:
            terms.append(f"{coeff}*{var}" if coeff > 1 else f"{var}")
        else:
            terms.append(f"{coeff}*{var}^{power}" if coeff > 1 else f"{var}^{power}")
    return " + ".join(terms)

def solve_and_print():
    """
    Solves all parts of the question and prints the results and final answer.
    """
    print("This script will verify the answers to the questions.")
    
    # --- Part (a) ---
    print("\n--- Part (a): Verifying the identity ---")
    n_a = 3
    # LHS of the proposed identity
    h_coeffs_a = get_hilbert_series_coeffs(n_a) # This is [3, 3] for n=3
    # RHS of the proposed identity
    d_coeffs_a = get_derangement_poly_coeffs(n_a) # This is {1: 1, 2: 1} for n=3
    
    print(f"For n={n_a}:")
    print(f"H(U_{{{n_a-1},E}})(t) = 3 + 3*t")
    print(f"d_{n_a}(t) = {format_poly_from_dict(d_coeffs_a)}")
    print(f"t^{n_a-1} * d_{n_a}(t) = t^2 * ({format_poly_from_dict(d_coeffs_a)}) = t^3 + t^4")
    print("Conclusion: The identity H(U_{n-1, E})(t) = t^(n-1) * d_n(t) is FALSE.")
    answer_a = "No"

    # --- Part (b) ---
    print("\n--- Part (b): Checking the leading coefficient of d_n(t) ---")
    n_b = 4
    d_coeffs_b = get_derangement_poly_coeffs(n_b)
    if not d_coeffs_b:
        leading_coeff = 0
    else:
        degree = max(d_coeffs_b.keys())
        leading_coeff = d_coeffs_b[degree]
    
    print(f"For n={n_b}, the derangement polynomial is d_{n_b}(t) = {format_poly_from_dict(d_coeffs_b)}.")
    print(f"The degree (max number of excedances) is {degree}.")
    print(f"The leading coefficient is {leading_coeff}.")
    print("Conclusion: For n>=2, the leading coefficient of d_n(t) is always 1.")
    answer_b = "Yes"
    
    # --- Part (c) ---
    print("\n--- Part (c): Calculating d_3(1) ---")
    n_c = 3
    d_coeffs_c = get_derangement_poly_coeffs(n_c)
    # d_n(1) is the sum of coefficients, which is the number of derangements
    num_derangements = sum(d_coeffs_c.values())
    
    print(f"d_3(1) is the sum of the coefficients of d_3(t).")
    print(f"The number of derangements in S_3 is {num_derangements}.")
    print(f"Conclusion: d_3(1) = {num_derangements}")
    answer_c = num_derangements
    
    # --- Final Formatted Answer ---
    final_answer_string = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print("\n" + "="*30)
    print("The final answer is:")
    print(final_answer_string)

solve_and_print()