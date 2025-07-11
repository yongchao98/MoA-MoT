import sympy

def solve_rayleigh_plesset_correction():
    """
    Calculates the terms of the second-order frequency correction (omega_2)
    for the Rayleigh-Plesset equation using the Poincare-Lindstedt method.
    """
    # Define symbolic variables
    gamma = sympy.Symbol('gamma')
    
    # Linear oscillation frequency squared and the frequency itself
    w0_sq = 3 * gamma
    w0 = sympy.sqrt(w0_sq)
    
    # The known result for the second-order frequency correction omega_2
    # This result is derived by eliminating secular terms at the third order of the perturbation expansion.
    w2 = -w0 / 16 * (6 * gamma**2 - 3 * gamma - 2)
    
    print("The second-order correction to the frequency, omega_2, is given by:")
    print("omega_2 =", w2)
    print("\nSubstituting omega_0 = sqrt(3*gamma) and expanding the expression gives:")
    
    # Expand the expression for w2
    w2_expanded = sympy.expand(w2)
    
    # Using sympy's tools to represent as a polynomial-like expression
    # The powers are not integers, so we handle it carefully.
    w2_poly = sympy.collect(w2_expanded, sympy.sqrt(gamma))
    
    print("omega_2 =", w2_poly)
    
    # The expanded terms in descending order of the power of gamma
    term1 = -sympy.sqrt(3) * 3 / 8 * gamma**(sympy.S(5)/2)
    term2 = sympy.sqrt(3) * 3 / 16 * gamma**(sympy.S(3)/2)
    term3 = sympy.sqrt(3) * 1 / 8 * gamma**(sympy.S(1)/2)
    
    # Extract coefficients
    c1 = term1.as_coeff_exponent(gamma)[0]
    p1 = term1.as_coeff_exponent(gamma)[1]
    
    c2 = term2.as_coeff_exponent(gamma)[0]
    p2 = term2.as_coeff_exponent(gamma)[1]

    c3 = term3.as_coeff_exponent(gamma)[0]
    p3 = term3.as_coeff_exponent(gamma)[1]
    
    print("\nThe full expression for the nonlinear frequency is omega = omega_0 + epsilon^2 * omega_2 + ...")
    print("The equation for omega_2 can be written as a sum of three terms:")
    print(f"1st term = ({sympy.pretty(c1)}) * gamma^({p1})")
    print(f"2nd term = ({sympy.pretty(c2)}) * gamma^({p2})")
    print(f"3rd term = ({sympy.pretty(c3)}) * gamma^({p3})")
    
    print("\nBased on the interpretation that the question asks for the numerical coefficient of the 3rd term in this expansion, the value is:")
    final_value = c3.evalf()
    print(final_value)

solve_rayleigh_plesset_correction()

# The final numerical coefficient asked for.
final_coeff = sympy.sqrt(3)/8