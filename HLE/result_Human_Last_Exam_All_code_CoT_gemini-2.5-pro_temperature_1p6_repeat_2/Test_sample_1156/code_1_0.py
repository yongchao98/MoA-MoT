import sympy

def solve():
    """
    This function finds and prints the normalised density of the invariant measure
    for the map T(x) = 1/x**(1/2) mod 1.
    """
    # The unnormalized density is of the form K*(2*sqrt(x) + 1)
    x = sympy.Symbol('x')
    unnormalized_rho = 2 * sympy.sqrt(x) + 1

    # To find the normalization constant K, we integrate from 0 to 1 and set the result to 1.
    # K * Integral(rho_unnormalized) = 1
    integral_val = sympy.integrate(unnormalized_rho, (x, 0, 1))
    
    # The normalization constant is the reciprocal of the integral's value.
    K = 1 / integral_val
    
    # The normalized density is K * unnormalized_rho
    normalized_rho = K * unnormalized_rho
    
    # Expressing the constant and function separately.
    coeff_val = K.evalf()
    term1_coeff = (K*2).evalf()
    
    print("The normalised density of the invariant measure is rho(x) = K * (2*sqrt(x) + 1)")
    print(f"The normalization constant K is 1/({integral_val}) = {K}")
    print(f"So, the density function is: rho(x) = {K} * (2*sqrt(x) + 1)")
    print(f"rho(x) = {K*2}*sqrt(x) + {K}")

solve()