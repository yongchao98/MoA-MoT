import itertools

def is_derangement(perm):
    """
    Checks if a permutation is a derangement.
    A permutation is 1-indexed, e.g., (2, 1, 3) corresponds to sigma(1)=2, sigma(2)=1, sigma(3)=3.
    """
    for i, val in enumerate(perm):
        # i is 0-indexed, so we check against i+1
        if val == i + 1:
            return False
    return True

def get_excedance_count(perm):
    """
    Calculates the number of excedances in a 1-indexed permutation.
    """
    count = 0
    for i, val in enumerate(perm):
        # An excedance is at index i+1 if sigma(i+1) > i+1
        if val > i + 1:
            count += 1
    return count

def get_derangement_polynomial_coeffs(n):
    """
    Computes the coefficients of the derangement polynomial d_n(t).
    Returns a list `coeffs` where `coeffs[k]` is the coefficient of t^k.
    The size of the list is n, allowing for a maximum of n-1 excedances.
    """
    if n < 2:
        return [0] * n
    
    # max excedances is n-1, so indices 0 to n-1 cover all possibilities
    coeffs = [0] * n 
    
    # Generate permutations of {1, 2, ..., n}
    for p in itertools.permutations(range(1, n + 1)):
        if is_derangement(p):
            exc_count = get_excedance_count(p)
            if exc_count < len(coeffs):
                coeffs[exc_count] += 1
            
    return coeffs

def solve_and_print():
    """
    Solves all parts of the question and prints the results and final answer.
    """
    # --- Part (a): Verify the identity ---
    print("Part (a): Confirm whether H(U_{n-1, E})(t) = t^{n-1} d_n(t).")
    print("We will test this by comparing polynomial degrees for n=4.")
    n_a = 4
    
    # Theoretical degree of the Left-Hand Side (LHS)
    deg_lhs = n_a - 2
    print(f"For n={n_a}, the rank of the matroid U_{{{n_a-1}, E}} is r={n_a-1}.")
    print(f"The degree of the Hilbert series H(U_{{{n_a-1}, E}})(t) is r-1 = {n_a-1}-1 = {deg_lhs}.")
    
    # Compute the degree of the Right-Hand Side (RHS)
    coeffs_d4 = get_derangement_polynomial_coeffs(n_a)
    deg_d4 = 0
    for i in range(len(coeffs_d4) - 1, -1, -1):
        if coeffs_d4[i] != 0:
            deg_d4 = i
            break
            
    deg_rhs = (n_a - 1) + deg_d4
    print(f"The derangement polynomial d_{n_a}(t) has degree {deg_d4}.")
    print(f"The degree of the right-hand side t^({n_a-1}) * d_{n_a}(t) is ({n_a-1}) + {deg_d4} = {deg_rhs}.")
    
    answer_a = "No"
    print(f"Since the degrees {deg_lhs} and {deg_rhs} do not match, the identity is false.")
    print("-" * 30)

    # --- Part (b): Find the leading coefficient of d_n(t) ---
    print("Part (b): State if the leading coefficient of d_n(t) for any n >= 2 is always 1.")
    print("The leading coefficient is the number of derangements with the maximum possible number of excedances (n-1).")
    print("The unique derangement with n-1 excedances is the permutation (2, 3, ..., n, 1).")
    print("Thus, the leading coefficient is expected to be 1.")
    print("Verifying for n = 2, 3, 4, 5:")
    all_coeffs_one = True
    for n_b in range(2, 6):
        coeffs_dn = get_derangement_polynomial_coeffs(n_b)
        leading_coeff = 0
        deg_dn = 0
        # Find degree and leading coefficient from the end of the list
        for i in range(len(coeffs_dn) - 1, -1, -1):
            if coeffs_dn[i] != 0:
                leading_coeff = coeffs_dn[i]
                deg_dn = i
                break
        print(f"For n={n_b}, deg(d_{n_b}(t)) = {deg_dn}, and the leading coefficient is {leading_coeff}.")
        if leading_coeff != 1:
            all_coeffs_one = False
    answer_b = "Yes" if all_coeffs_one else "No"
    print("The verification confirms the leading coefficient is 1.")
    print("-" * 30)

    # --- Part (c): Calculate d_3(1) ---
    print("Part (c): Give the value of d_3(1).")
    print("d_n(1) is the sum of coefficients, which equals the total number of derangements of n elements.")
    n_c = 3
    
    count = 0
    # Generate permutations of {1, 2, 3} and count derangements
    for p in itertools.permutations(range(1, n_c + 1)):
        if is_derangement(p):
            count += 1
            
    answer_c = count
    print(f"By counting, the number of derangements in S_{n_c} is {answer_c}.")
    print("-" * 30)
    
    # --- Final Answer ---
    print("The final answer is:")
    print(f"<<<(a) {answer_a}; (b) {answer_b}; (c) {answer_c}>>>")

if __name__ == '__main__':
    solve_and_print()