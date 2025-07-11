import sympy

def solve_vector_sum_limit():
    """
    This function symbolically calculates the limit of n*P(n) as described in the problem.
    """
    # Define the symbol n for the total number of vectors.
    n = sympy.Symbol('n', positive=True, real=True)

    # Step 1: Calculate the elements of the covariance matrix of the sum S = (Sx, Sy).
    # We have n/3 vectors of each type: v_A=(1,0), v_B=(0.5, sqrt(3)/2), v_C=(-0.5, sqrt(3)/2).
    # Var(Sx) = sum of v_ix^2 over all n vectors.
    var_sx = (n/3) * (1**2) + (n/3) * (sympy.S(1)/2)**2 + (n/3) * (sympy.S(-1)/2)**2
    # Var(Sy) = sum of v_iy^2 over all n vectors.
    var_sy = (n/3) * (0**2) + (n/3) * (sympy.sqrt(3)/2)**2 + (n/3) * (sympy.sqrt(3)/2)**2
    # Cov(Sx, Sy) = sum of v_ix * v_iy over all n vectors.
    cov_sxy = (n/3) * (1 * 0) + (n/3) * (sympy.S(1)/2 * sympy.sqrt(3)/2) + (n/3) * (sympy.S(-1)/2 * sympy.sqrt(3)/2)

    # Simplify these expressions.
    var_sx_simplified = sympy.simplify(var_sx)
    var_sy_simplified = sympy.simplify(var_sy)
    cov_sxy_simplified = sympy.simplify(cov_sxy)

    print("Step 1: Determine the parameters of the approximating normal distribution.")
    print(f"The mean of the sum S is (0, 0).")
    print(f"The variance Var(Sx) is: {var_sx_simplified}")
    print(f"The variance Var(Sy) is: {var_sy_simplified}")
    print(f"The covariance Cov(Sx, Sy) is: {cov_sxy_simplified}")
    print("-" * 40)

    # Step 2: Calculate the probability P(n).
    # P(n) is the integral of the normal PDF over the disk x^2+y^2 <= 2.
    # In polar coordinates, this integral simplifies to the one below.
    r = sympy.Symbol('r', positive=True)
    # The integral is derived from integrating the PDF f(r) = (2*pi*r) * (1/(pi*n)) * exp(-r^2/n)
    integrand = (2 / n) * r * sympy.exp(-r**2 / n)
    
    # Perform the definite integration from r=0 to r=sqrt(2).
    Pn_expression = sympy.integrate(integrand, (r, 0, sympy.sqrt(2)))
    
    print("Step 2: Calculate the probability P(n).")
    print(f"The probability P(n) for large n is given by the expression: {Pn_expression}")
    print("-" * 40)

    # Step 3: Calculate the final limit of n * P(n).
    expression_to_limit = n * Pn_expression
    final_limit = sympy.limit(expression_to_limit, n, sympy.oo)

    print("Step 3: Compute the limit of n * P(n) as n -> infinity.")
    print(f"We want to compute the limit of: {expression_to_limit}")
    print(f"The final equation is limit(n * P(n)) = L, where L is the result.")
    # Here we output the numbers and operators in the final equation as requested.
    print(f"limit(n * ({Pn_expression})) as n->oo = {final_limit}")


solve_vector_sum_limit()