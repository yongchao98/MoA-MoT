import sympy

def solve_part_a_counterexample():
    """
    This function provides a computational verification of the counterexample for part (a).
    It sets up the specific family of sets and polynomials described in the analysis,
    and then checks for linear independence using symbolic computation.
    """
    # Setup for the counterexample: n=3, s=2
    n = 3
    L = {0, 1}
    # The family F is all 2-element subsets of {1, 2, 3}
    F = [{1, 2}, {1, 3}, {2, 3}]
    m = len(F)

    # Symbolic variables
    x = sympy.symbols(f'x_1:{n+1}')
    c = sympy.symbols(f'c_1:{m+1}')

    print("--- Counterexample Verification for Part (a) ---")
    print(f"Let n = {n}, L = {L}. This gives s = |L| = {len(L)}.")
    print(f"The condition s > floor(n/2) becomes {len(L)} > {n//2}, which is True.")
    print(f"Let F be the family of all 2-element subsets of [n]: {F}")
    print("-" * 50)

    # Construct characteristic vectors
    vectors = []
    for s_i in F:
        v = [0] * n
        for i in s_i:
            v[i-1] = 1
        vectors.append(v)

    # Construct polynomials P_i(x)
    polys = []
    print("The constructed polynomials P_i(x) are:")
    for i, v_i in enumerate(vectors):
        F_i = F[i]
        dot_product = sum(x[j] * v_i[j] for j in range(n))
        
        # P_i(x) = prod_{l in L, l < |F_i|} (<x,v_i> - l)
        # Here |F_i|=2 for all i, so the condition l < |F_i| holds for all l in L.
        p = 1
        for l_k in L:
            p *= (dot_product - l_k)
        
        expanded_p = p.expand()
        polys.append(expanded_p)
        print(f"P_{i+1}(x) = {sympy.pretty(expanded_p, use_unicode=True)}\n")

    # Check for linear independence
    # Set up the equation: c_1*P_1 + c_2*P_2 + ... = 0
    linear_combination = sum(c[i] * polys[i] for i in range(m))

    print("-" * 50)
    print("To check for linear independence, we solve the equation c_1*P_1 + c_2*P_2 + c_3*P_3 = 0.")
    # The equation must hold for all x, so the coefficient of each monomial must be zero.
    
    # Collect all unique monomials from the expanded sum
    all_monomials = set()
    for p in polys:
        all_monomials.update(p.as_poly(*x).monoms())
    
    # Generate system of equations for coefficients c_i
    equations = []
    for mono_tuple in all_monomials:
        mono = x[0]**mono_tuple[0] * x[1]**mono_tuple[1] * x[2]**mono_tuple[2]
        coeff = linear_combination.coeff(mono)
        if coeff != 0:
            equations.append(coeff)
    
    print("This yields the following system of linear equations for c_1, c_2, c_3:")
    for eq in equations:
        print(f"  {eq} = 0")
        
    # Solve the system
    solution = sympy.solve(equations, c)

    print("\nSolving this system gives:")
    print(f"  {solution}")

    # Check if the only solution is the trivial one
    is_trivial = all(value == 0 for value in solution.values())
    if is_trivial:
        print("\nConclusion: The only solution is c_1=c_2=c_3=0. The polynomials are linearly independent.")
        print("This provides a counterexample to the statement in (a).")
    else:
        print("\nConclusion: A non-trivial solution exists. The polynomials are linearly dependent.")

if __name__ == '__main__':
    solve_part_a_counterexample()