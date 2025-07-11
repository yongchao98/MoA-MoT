import sympy

def check_linear_independence():
    """
    This function constructs a counterexample for question (a) and checks for linear independence.
    We choose n=3, L={0, 1}. This gives s=2, which satisfies s > floor(n/2) (2 > 1).
    The family F is an ordered L-intersecting family of subsets of [3] (using 1-based indexing for sets).
    F = {{1,3}, {2,3}, {1,2}}.
    It is ordered because we can set r=2. 
    - n=3 is in F1={1,3} and F2={2,3}.
    - n=3 is not in F3={1,2}.
    - |F1|=|F2|=|F3|=2, so the size-ordering condition |Fi|<=|Fj| for i<j holds.
    It is L-intersecting because the pairwise intersections are |F1_cap_F2|=1, |F1_cap_F3|=1, |F2_cap_F3|=1, all in L.

    The polynomials Pi are defined as Prod_{k: l_k < |Fi|} ( <x, vi> - l_k ).
    Since |Fi|=2 for all i, and L={0, 1}, the condition l_k < 2 includes both elements of L.
    So, Pi(x) = (<x, vi> - 0) * (<x, vi> - 1).
    """

    # Define symbolic variables for x = (x1, x2, x3)
    x1, x2, x3 = sympy.symbols('x1 x2 x3')
    
    # Define the sets F in terms of 1-based indices
    F = [{1, 3}, {2, 3}, {1, 2}]
    n = 3

    # Define characteristic vectors vi for each set Fi in F
    # v1 for {1,3} -> (1, 0, 1)
    # v2 for {2,3} -> (0, 1, 1)
    # v3 for {1,2} -> (1, 1, 0)
    v = [(1, 0, 1), (0, 1, 1), (1, 1, 0)]
    x_vec = [x1, x2, x3]

    polynomials = []
    print("Constructing polynomials for the counterexample (n=3, s=2, L={0,1}):")
    # L = {0, 1}, so l1=0, l2=1
    # For all sets Fi, |Fi|=2, so we multiply over l_k < 2, which means l_k=0 and l_k=1.
    for i in range(len(F)):
        # Calculate the dot product <x, vi>
        dot_product = sum(x_j * v_ij for x_j, v_ij in zip(x_vec, v[i]))
        
        # Construct the polynomial Pi = (<x, vi> - 0) * (<x, vi> - 1)
        p = dot_product * (dot_product - 1)
        polynomials.append(p)
        print(f"P{i+1}(x1,x2,x3) = {sympy.expand(p)}")

    # Check for linear independence.
    # We set up the equation c1*P1 + c2*P2 + c3*P3 = 0 and solve for c1, c2, c3.
    # The polynomials are linearly independent if the only solution is c1=c2=c3=0.
    
    # Create symbols for coefficients
    c1, c2, c3 = sympy.symbols('c1 c2 c3')
    
    # Form the linear combination
    linear_combination = c1 * polynomials[0] + c2 * polynomials[1] + c3 * polynomials[2]
    
    # Expand the combination
    expanded_lc = sympy.expand(linear_combination)
    
    # Collect coefficients for each monomial in x1, x2, x3
    # The linear combination must be the zero polynomial, so all coefficients must be zero.
    coeffs = sympy.Poly(expanded_lc, c1, c2, c3).coeffs()
    
    # The system of equations is derived from setting the coefficient of each x-monomial to zero
    # For example, coefficient of x1**2 must be zero, etc.
    eqs = [sympy.Eq(expanded_lc.coeff(m), 0) for m in expanded_lc.free_symbols if m not in {c1, c2, c3}]

    # A simpler way: get all coefficients of the expanded polynomial with respect to x variables
    poly_in_x = sympy.Poly(expanded_lc, x1, x2, x3)
    eq_system = poly_in_x.coeffs()

    # Solve the system of linear equations for c1, c2, c3
    solution = sympy.solve(eq_system, (c1, c2, c3), dict=True)

    print("\nSolving the equation c1*P1 + c2*P2 + c3*P3 = 0 for coefficients c1, c2, c3.")
    # The only solution is the trivial one {c1: 0, c2: 0, c3: 0}
    if solution and all(val == 0 for val in solution[0].values()):
        print("The only solution is c1=0, c2=0, c3=0.")
        print("Result: The polynomials are linearly independent.")
    else:
        print("A non-trivial solution exists.")
        print("Result: The polynomials are linearly dependent.")

    print("\n---\nFinal Answer to Question (a):")
    print("The statement is FALSE. The provided counterexample shows a case where the polynomials are linearly independent.")


if __name__ == '__main__':
    check_linear_independence()
