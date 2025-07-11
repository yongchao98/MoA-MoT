import math

def solve_integral_approximation():
    """
    This function calculates and displays the analytical approximation for the given integral
    in the small epsilon regime.
    """
    
    # The integral is I(eps) = integral from 0 to 15 of 1/(eps + g(x)) dx
    # where g(x) = 9.0*x^5.0 + 5.0*x^6.0 + 9.0*x^8.0.
    
    # For small epsilon, the integral is dominated by the behavior near x=0.
    # We approximate g(x) by its leading term for small x, which is G * x^P.
    G = 9.0
    P = 5.0
    
    print("### Derivation of the Analytical Formula ###\n")
    print(f"The integral is approximated by I(eps) ~= Integral(0 to inf) [1 / (eps + {G}*x^{P})] dx.")
    print("This integral has an analytical solution of the form: C * eps^p\n")

    # Calculate the exponent 'p'
    # p = 1/P - 1
    p_exponent = 1.0/P - 1.0
    
    print(f"The exponent 'p' is calculated as: p = 1/P - 1 = 1/{P} - 1 = {p_exponent}\n")

    # Calculate the coefficient 'C'
    # C = (pi / (P * sin(pi/P))) * G^(-1/P)
    # This can be rewritten for clarity as: pi / (P * G^(1/P) * sin(pi/P))
    pi_val = math.pi
    sin_pi_over_P = math.sin(pi_val / P)
    G_pow_1_over_P = math.pow(G, 1.0/P)
    
    C_coefficient = pi_val / (P * G_pow_1_over_P * sin_pi_over_P)

    print("The coefficient 'C' is calculated using the formula: C = pi / (P * G^(1/P) * sin(pi/P))")
    print("Plugging in the numbers:")
    print(f"pi = {pi_val}")
    print(f"P = {P}")
    print(f"G = {G}")
    print(f"sin(pi/P) = sin(pi/{P}) = {sin_pi_over_P}")
    print(f"G^(1/P) = {G}^(1/{P}) = {G_pow_1_over_P}\n")
    
    print("### Final Analytical Formula ###\n")
    print("The final approximate formula, with all numerical values plugged in, is:")
    print(f"I(epsilon) approx ({pi_val:.4f} / ({P} * {G_pow_1_over_P:.4f} * {sin_pi_over_P:.4f})) * epsilon^({p_exponent})")
    print("\nWhich simplifies to:")
    print(f"I(epsilon) approx {C_coefficient} * epsilon^({p_exponent})")

solve_integral_approximation()