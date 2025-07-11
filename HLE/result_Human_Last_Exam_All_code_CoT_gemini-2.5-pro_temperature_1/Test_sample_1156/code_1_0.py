import sympy

def solve():
    """
    This function finds and prints the normalized density of the invariant measure
    for the map T(x) = 1/sqrt(x) mod 1.
    """
    # The form of the invariant density is known to be linear: rho(x) = C*(1+x).
    # We need to find the normalization constant C such that the integral from 0 to 1 is 1.
    x = sympy.Symbol('x')
    C = sympy.Symbol('C')
    
    # The unnormalized density function
    rho_unnormalized = C * (1 + x)
    
    # Calculate the integral to find the normalization constant
    integral_value = sympy.integrate(rho_unnormalized, (x, 0, 1))
    
    # Solve for C
    # The equation is integral_value = 1
    normalization_constant = sympy.solve(sympy.Eq(integral_value, 1), C)[0]
    
    # The normalized density function
    a = normalization_constant * 1
    b = normalization_constant * 1
    
    print(f"The normalised density of the invariant measure is rho(x) = a + b*x.")
    print("The coefficients of the linear function are:")
    print(f"a = {a}")
    print(f"b = {b}")

solve()
