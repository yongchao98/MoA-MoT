import sympy

def solve_polynomial_extremum():
    """
    Solves for the infimum and supremum of the given expression involving polynomial roots and coefficients.
    """
    # Set up the symbolic variables
    x = sympy.Symbol('x')
    f = x**3 + x**2 + x - 1

    # --- Part 1: Asymptotic analysis ---
    # Find roots
    roots_map = sympy.roots(f, multiple=True)
    roots = list(roots_map.keys())
    
    xi = None
    complex_roots = []
    for r in roots:
        if r.is_real:
            xi = r
        else:
            complex_roots.append(r)

    # Ensure xi1 has a positive imaginary part
    if complex_roots[0].as_real_imag()[1] > 0:
        xi1 = complex_roots[0]
    else:
        xi1 = complex_roots[1]

    # Calculate the vector K for asymptotic analysis
    f_prime = sympy.diff(f, x)
    K_vec = sympy.Matrix([1/xi1, 1 + xi1, 1]) / f_prime.subs(x, xi1)

    # Calculate the two key dot products
    K_dot_K = K_vec.dot(K_vec)
    K_H_sq = K_vec.dot(K_vec.conjugate())

    # Asymptotic infimum and supremum from theoretical derivation
    inf_limit = sympy.re(2 * (K_H_sq - sympy.Abs(K_dot_K)))
    sup_limit = sympy.re(2 * (K_H_sq + sympy.Abs(K_dot_K)))

    # --- Part 2: Calculation for small n ---
    # Initialize coefficients for P_1(x) = x
    a, b, c = 0, 1, 0
    xi_val = xi.evalf(100)
    
    min_val = sup_limit.evalf(100) # Initialize min with a large number
    max_val = inf_limit.evalf(100) # Initialize max with a small number

    num_iterations = 50
    for n in range(1, num_iterations + 1):
        # Calculate E_n = |xi^n| * (a^2 + b^2 + c^2)
        # For real 0 < xi < 1, |xi^n| = xi^n
        current_val = (xi_val**n) * (a**2 + b**2 + c**2)
        
        if n == 1:
            min_val = current_val
            max_val = current_val
        else:
            if current_val < min_val:
                min_val = current_val
            if current_val > max_val:
                max_val = current_val
            
    # Update coefficients for P_{n+1}(x) using the recurrence relation
    # a_{n+1} = c_n, b_{n+1} = a_n - c_n, c_{n+1} = b_n - c_n
    # Here, (a,b,c) are coefficients for P_n(x)
        a_next = c
        b_next = a - c
        c_next = b - c
        a, b, c = a_next, b_next, c_next

    # The final infimum and supremum are the min/max of the computed initial values
    # and the asymptotic limits, as the sequence converges quickly.
    final_inf = min(min_val, inf_limit.evalf(100))
    final_sup = max(max_val, sup_limit.evalf(100))

    # Output the final answer per instructions
    print(f"inf_n |P_n(xi)|(a_n^2+b_n^2+c_n^2) = {final_inf}")
    print(f"sup_n |P_n(xi)|(a_n^2+b_n^2+c_n^2) = {final_sup}")

solve_polynomial_extremum()