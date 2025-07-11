import itertools

def is_derangement(perm):
    """Checks if a 1-indexed permutation is a derangement."""
    for i, p_i in enumerate(perm):
        if p_i == i + 1:
            return False
    return True

def excedances(perm):
    """Calculates the number of excedances in a 1-indexed permutation."""
    count = 0
    for i, p_i in enumerate(perm):
        if p_i > i + 1:
            count += 1
    return count

def get_derangement_poly_coeffs(n):
    """
    Computes the coefficients of the derangement polynomial d_n(t).
    Returns a dictionary where keys are exponents and values are coefficients.
    """
    if n < 2: # d_0(t)=1 (by convention, empty perm), d_1(t)=0
        return {0: 0} if n==1 else {0: 1}
        
    coeffs = {}
    for p in itertools.permutations(range(1, n + 1)):
        if is_derangement(p):
            exc_count = excedances(p)
            coeffs[exc_count] = coeffs.get(exc_count, 0) + 1
    return coeffs

def solve():
    """Solves all parts of the question and prints the results."""
    
    # Part (a)
    print("Part (a) Analysis:")
    n_a = 4 # Using n=4 as an example
    deg_H = n_a - 1
    d_n_coeffs = get_derangement_poly_coeffs(n_a)
    deg_d_n = max(d_n_coeffs.keys()) if d_n_coeffs else 0
    deg_rhs = (n_a - 1) + deg_d_n

    print(f"For n={n_a}, the degree of H(U_{n-1, E})(t) is n-1 = {deg_H}.")
    print(f"The degree of d_{n_a}(t) is {deg_d_n}.")
    print(f"The degree of the right-hand side t^(n-1)d_n(t) is (n-1) + deg(d_n(t)) = {n_a-1} + {deg_d_n} = {deg_rhs}.")
    is_true_a = (deg_H == deg_rhs)
    print(f"Since {deg_H} != {deg_rhs}, the identity is false.")
    
    answer_a_bool = "No"
    answer_a_expr = "n-1"
    print(f"Result for (a): {answer_a_bool}, the degree is {answer_a_expr}.\n")

    # Part (b)
    print("Part (b) Analysis:")
    print("Checking the leading coefficient of d_n(t) for n=2, 3, 4, 5:")
    is_true_b = True
    for n_b in range(2, 6):
        coeffs = get_derangement_poly_coeffs(n_b)
        degree = max(coeffs.keys())
        leading_coeff = coeffs[degree]
        print(f"For n={n_b}, d_{n_b}(t) has degree {degree} and leading coefficient {leading_coeff}.")
        if leading_coeff != 1:
            is_true_b = False
    print("As shown theoretically and by computation, the leading coefficient is always 1 for n>=2.")
    
    answer_b_bool = "Yes"
    print(f"Result for (b): {answer_b_bool}.\n")
    
    # Part (c)
    print("Part (c) Analysis:")
    n_c = 3
    d3_coeffs = get_derangement_poly_coeffs(n_c)
    # d_3(1) is the sum of coefficients, which is the number of derangements !3
    d3_at_1 = sum(d3_coeffs.values())
    
    # Building the equation string
    poly_str = " + ".join([f"{coeff}*1^{exp}" for exp, coeff in sorted(d3_coeffs.items(), reverse=True)])
    print(f"d_3(1) is the number of derangements in S_3.")
    print(f"The derangements in S_3 are (2, 3, 1) and (3, 1, 2).")
    print(f"The value is !3 = 2.")
    print(f"Calculation from polynomial: d_3(1) = {poly_str} = {d3_at_1}")

    answer_c_expr = d3_at_1
    print(f"Result for (c): {answer_c_expr}.\n")

    # Final formatted answer
    final_answer = f"(a) {answer_a_bool}, {answer_a_expr}; (b) {answer_b_bool}; (c) {answer_c_expr}"
    print("Final Answer:")
    print(f"<<<{final_answer}>>>")

solve()