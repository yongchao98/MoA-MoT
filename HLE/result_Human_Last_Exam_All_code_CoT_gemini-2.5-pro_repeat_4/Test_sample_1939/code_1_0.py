import sympy

def solve():
    """
    Computes the Ehrhart polynomial and its roots for the given polytope
    for dimensions d=2, 3, and 4.
    """
    z = sympy.Symbol('z')

    for d in range(2, 5):
        print(f"--- Analyzing for dimension d={d} ---")
        n = sympy.Symbol('n')
        k = sympy.Symbol('k')

        # Formula for the number of integer points in the k-th slice
        # N_n(k) = sum_{j=0}^{d-1} C(d-1, j) * C(n-k+j, j) * C(k, d-1-j)
        j = sympy.Symbol('j')
        term = sympy.binomial(d-1, j) * sympy.binomial(n - k + j, j) * sympy.binomial(k, d - 1 - j)
        N_k = sympy.summation(term, (j, 0, d-1))

        # The Ehrhart polynomial p(n) is the sum of N_k from k=0 to n
        # We replace n with z for the polynomial form
        p_poly_expr = sympy.summation(N_k.doit(), (k, 0, n))
        p_poly_expr = p_poly_expr.expand()
        p_poly = sympy.Poly(p_poly_expr, n)
        
        # Substitute n with z for the final polynomial
        p_z = p_poly.subs(n, z)
        
        print(f"The Ehrhart polynomial p(z) is: {p_z.as_expr()}")

        # Find the roots of the polynomial
        try:
            roots = sympy.roots(p_z)
            print("The roots of p(z) are:")
            # We print the roots in a more readable format
            root_list = []
            for r, multiplicity in roots.items():
                for _ in range(multiplicity):
                    root_list.append(r)
            
            for r in sorted(root_list, key=lambda x: (sympy.re(x), sympy.im(x))):
                print(f"  {r}")

            # Check the properties based on the roots
            real_parts = [sympy.re(r) for r in root_list]
            is_real_part_minus_1 = all(rp == -1 for rp in real_parts)
            
            print(f"\nConclusion for d={d}:")
            if is_real_part_minus_1:
                print("Every root of p(z) has real part -1.")
            else:
                print("Not every root of p(z) has real part -1.")

        except Exception as e:
            print(f"Could not compute roots for d={d}: {e}")
        print("-" * 30)

solve()