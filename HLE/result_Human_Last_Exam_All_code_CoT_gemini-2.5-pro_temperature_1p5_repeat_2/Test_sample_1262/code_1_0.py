import itertools

def get_derangement_poly_coeffs(n):
    """
    Calculates the coefficients of the n-th derangement polynomial d_n(t).

    The definition of permutation, derangement, and excedance are 1-indexed.
    Python lists are 0-indexed, so we use elements {0, ..., n-1}.
    A permutation p of {0, ..., n-1} is a derangement if p[i] != i.
    An excedance for p is an index i such that p[i] > i.
    The definitions in the problem are for {1, ..., n}.
    sigma(i) > i is equivalent to p[i-1] > i-1.
    
    Args:
        n (int): The number of elements.

    Returns:
        list: A list of coefficients of d_n(t), where coeffs[k] is the
              coefficient of t^k.
    """
    if n == 0:
        return [1]
    if n == 1:
        return [0]

    # The maximum number of excedances is n-1. The polynomial degree is at most n-1.
    # coeffs[k] will store the number of derangements with k excedances.
    coeffs = [0] * n 
    
    # Generate all permutations of {0, 1, ..., n-1}
    elements = range(n)
    for p in itertools.permutations(elements):
        # Check if it is a derangement
        is_derangement = all(p[i] != i for i in range(n))
        
        if is_derangement:
            # Count excedances
            excedances = sum(1 for i in range(n) if p[i] > i)
            if excedances < n:
                coeffs[excedances] += 1
            
    return coeffs

def get_poly_degree(coeffs):
    """Finds the degree of a polynomial given its coefficients."""
    for i in range(len(coeffs) - 1, -1, -1):
        if coeffs[i] != 0:
            return i
    return -1 # Degree of zero polynomial

def main():
    """
    Solves the three parts of the problem.
    """
    print("--- Analysis Start ---")

    # Part (a)
    n_a = 4
    deg_H = n_a - 1
    d_n_coeffs_a = get_derangement_poly_coeffs(n_a)
    deg_d_n = get_poly_degree(d_n_coeffs_a)
    deg_RHS = (n_a - 1) + deg_d_n
    
    print(f"\n(a) Analysis for n = {n_a}:")
    print(f"The degree of H(U_{n-1,E})(t) is the rank of the matroid, which is n-1 = {deg_H}.")
    print(f"The computed degree of the derangement polynomial d_{n_a}(t) is {deg_d_n}.")
    print(f"The degree of the expression t^(n-1) * d_n(t) is (n-1) + deg(d_n(t)) = {n_a-1} + {deg_d_n} = {deg_RHS}.")
    print(f"Since the degrees ({deg_H} and {deg_RHS}) are not equal for n={n_a}, the identity is False.")

    # Part (b)
    print("\n(b) Analysis for the leading coefficient of d_n(t):")
    is_always_one = True
    for n_b in range(2, 9):
        coeffs_b = get_derangement_poly_coeffs(n_b)
        # Degree is n-1, so leading coefficient is at index n-1
        leading_coeff = coeffs_b[n_b - 1]
        print(f"For n = {n_b}, the leading coefficient of d_n(t) (coefficient of t^{n_b-1}) is {leading_coeff}.")
        if leading_coeff != 1:
            is_always_one = False
    if is_always_one:
        print("The leading coefficient is 1 for all tested n >= 2. This supports the 'Yes' answer.")
    else:
        print("The leading coefficient is not always 1.")
        
    # Part (c)
    n_c = 3
    coeffs_c = get_derangement_poly_coeffs(n_c)
    val_d3_1 = sum(coeffs_c)
    
    poly_terms = []
    for i, coeff in enumerate(coeffs_c):
        if coeff > 0:
            poly_terms.append(f"{coeff}*t^{i}")
    poly_str = " + ".join(poly_terms)

    coeff_values = [str(c) for c in coeffs_c if c > 0]
    
    print(f"\n(c) Calculation of d_3(1):")
    print(f"The derangement polynomial is d_3(t) = {poly_str}.")
    print(f"Evaluating at t=1 gives d_3(1), which is the sum of the coefficients.")
    print(f"The value is the total number of derangements for n=3.")
    print(f"d_3(1) = {' + '.join(coeff_values)} = {val_d3_1}")
    
    print("\n--- Analysis End ---")

if __name__ == "__main__":
    main()
