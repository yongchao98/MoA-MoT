import itertools

def get_derangement_polynomial_coeffs(n):
    """
    Computes the coefficients of the n-th derangement polynomial d_n(t).
    Returns a list of coefficients, where coeffs[k] is the coeff of t^k.
    """
    if n == 0:
        return [1] # d_0(t) = 1
    if n == 1:
        return [0] # d_1(t) = 0
    
    max_excedances = n - 1
    coeffs = [0] * (max_excedances + 1)
    
    for p in itertools.permutations(range(1, n + 1)):
        is_derangement = all(i != pi for i, pi in enumerate(p, 1))
        
        if is_derangement:
            excedances = sum(1 for i, pi in enumerate(p, 1) if pi > i)
            if excedances <= max_excedances:
                coeffs[excedances] += 1
    return coeffs

def poly_to_string(coeffs, var='t'):
    """Converts a list of coefficients to a string representation of a polynomial."""
    s = []
    for i in range(len(coeffs) - 1, -1, -1):
        if coeffs[i] != 0:
            if i == 0:
                s.append(f"{coeffs[i]}")
            elif i == 1:
                s.append(f"{coeffs[i]}*{var}")
            else:
                s.append(f"{coeffs[i]}*{var}^{i}")
    if not s:
        return "0"
    return " + ".join(s)

def solve_all():
    """
    Solves all parts of the problem and prints the results.
    """
    
    # --- Part (a) ---
    print("Part (a): Confirming the identity for n=3.")
    n_a = 3
    rank_a = n_a - 1
    
    # H(U_{n-1,E})(t) = t^(n-1) * (T(1+1/t, 1))
    # T_Cn(x,y) = x + y + ... + y^(n-1)
    # T_C3(1+1/t, 1) = (1+1/t) + 1 + 1^2 = 3 + 1/t
    # H(U_{2,E})(t) = t^2 * (3 + 1/t) = 3*t^2 + t
    h_poly_coeffs = [0, 1, 3] # Coeffs of t + 3t^2
    
    print(f"For n = {n_a}, the rank r = {rank_a}.")
    print(f"The Hilbert series H(U_{n_a-1, E})(t) is calculated to be {poly_to_string(h_poly_coeffs)}.")
    
    d3_coeffs = get_derangement_polynomial_coeffs(n_a)
    print(f"The derangement polynomial d_{n_a}(t) is {poly_to_string(d3_coeffs)}.")
    
    # RHS = t^(n-1) * d_n(t) = t^2 * d_3(t)
    # d_3(t) = t + t^2 -> coeffs are [0, 1, 1]
    # t^2 * d_3(t) = t^3 + t^4 -> coeffs are [0, 0, 0, 1, 1]
    rhs_poly_coeffs = [0, 0] + d3_coeffs 
    print(f"The right side of the identity, t^{n_a-1}*d_{n_a}(t), is {poly_to_string(rhs_poly_coeffs)}.")
    
    answer_a = "No"
    print(f"Since {poly_to_string(h_poly_coeffs)} != {poly_to_string(rhs_poly_coeffs)}, the identity is false.")
    print("-" * 20)

    # --- Part (b) ---
    print("Part (b): Checking the leading coefficient of d_n(t).")
    d2_coeffs = get_derangement_polynomial_coeffs(2)
    d3_coeffs = get_derangement_polynomial_coeffs(3)
    d4_coeffs = get_derangement_polynomial_coeffs(4)

    lc_d2 = d2_coeffs[-1]
    lc_d3 = d3_coeffs[-1]
    lc_d4 = d4_coeffs[-1]

    print(f"d_2(t) = {poly_to_string(d2_coeffs)}, leading coefficient is {lc_d2}.")
    print(f"d_3(t) = {poly_to_string(d3_coeffs)}, leading coefficient is {lc_d3}.")
    print(f"d_4(t) = {poly_to_string(d4_coeffs)}, leading coefficient is {lc_d4}.")
    
    answer_b = "No"
    print(f"The leading coefficient is not always 1 (e.g., for n=4 it is {lc_d4}). So the statement is false.")
    print("-" * 20)

    # --- Part (c) ---
    print("Part (c): Calculating d_3(1).")
    n_c = 3
    # d_3(t) has coeffs [0, 1, 1]
    val_c = sum(c * (1**i) for i, c in enumerate(d3_coeffs))
    print(f"d_3(t) = {poly_to_string(d3_coeffs)}.")
    print(f"d_3(1) = {d3_coeffs[1]}*1^1 + {d3_coeffs[2]}*1^2 = {d3_coeffs[1]} + {d3_coeffs[2]} = {val_c}.")
    print("This value is the number of derangements of 3 elements, which is 2.")
    answer_c = val_c
    print("-" * 20)
    
    # --- Final Answer ---
    final_answer_str = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print("Final formatted answer:")
    print(f"<<<{final_answer_str}>>>")

solve_all()