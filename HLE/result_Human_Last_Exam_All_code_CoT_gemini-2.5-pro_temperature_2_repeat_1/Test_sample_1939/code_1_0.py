import sympy

def solve_poly_properties():
    """
    This function analyzes the Ehrhart polynomial for the given polytope
    for dimensions d = 2, 3, and 4, and checks the given properties.
    """
    z = sympy.Symbol('z')

    for d in range(2, 5):
        print(f"--- Analysis for d = {d} ---")

        # The Ehrhart polynomial is p(z) = (z+1) * binomial(z+d-1, d-1)
        p = (z + 1) * sympy.binomial(z + d - 1, d - 1)
        p_expanded = sympy.expand(p)

        print(f"The Ehrhart polynomial p(z) is: {p_expanded}")

        # Find the roots of the polynomial
        # sympy.roots returns a dictionary of {root: multiplicity}
        roots_dict = sympy.roots(p, z)
        roots = list(roots_dict.keys())
        
        # sympy roots can be non-python types, so convert them to check properties
        real_parts = [sympy.re(r) for r in roots]
        is_real = [sympy.im(r) == 0 for r in roots]
        
        print(f"The roots of the polynomial are: {roots_dict}")

        # Check the answer choices
        # A. Every root of p has real part -1.
        choice_A = all(rp == -1 for rp in real_parts)
        print(f"A. Every root has real part -1: {choice_A}")

        # B. Every root of p is real.
        choice_B = all(is_real)
        print(f"B. Every root is real: {choice_B}")

        # C. The coefficients of p sum exactly d.
        # Sum of coefficients is p(1)
        sum_of_coeffs = p.subs(z, 1)
        choice_C = (sum_of_coeffs == d)
        print(f"C. The sum of coefficients is {sum_of_coeffs}, which should be equal to d={d}: {choice_C}")
        
        # D. The coefficients of p sum exactly d choose d/2.
        # For non-integer d/2, this is not well-defined.
        # We can test for even d.
        if d % 2 == 0:
            d_choose_d_half = sympy.binomial(d, d // 2)
            choice_D = (sum_of_coeffs == d_choose_d_half)
            print(f"D. The sum of coefficients is {sum_of_coeffs}, which should be equal to C(d,d/2)={d_choose_d_half}: {choice_D}")
        else:
            print(f"D. C(d, d/2) is not well-defined for odd d={d}. Statement is likely false.")
            
        # E. Every root of p has real part -1/2.
        choice_E = all(rp == -1/2 for rp in real_parts)
        print(f"E. Every root has real part -1/2: {choice_E}")
        print("-" * (20 + len(str(d))))

if __name__ == '__main__':
    solve_poly_properties()
