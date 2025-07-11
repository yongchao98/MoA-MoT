import sympy

def solve_polymer_force():
    """
    This function derives and prints the force law for a thermally isolated,
    freely jointed chain polymer.

    The derivation follows these steps:
    1.  The system is thermally isolated and stretched slowly, so its total entropy S is constant (dS = 0).
    2.  Total entropy S = S_configurational + S_kinetic.
        S_conf = C1 - (3 * k_B * x**2) / (2 * n * l**2)
        S_kin = C2 + (3*n/2) * k_B * log(K)
    3.  Setting dS = 0 yields a differential equation relating kinetic energy K and extension x:
        dK/K = (2 * x / (n**2 * l**2)) dx
    4.  The force F is the work done per unit extension, F = dK/dx.
        This gives F = (2 * K * x) / (n**2 * l**2).
    5.  Solving the differential equation from step 3 for K, with initial condition K(0) = E(0), gives:
        K(x) = E(0) * exp(x**2 / (n**2 * l**2))
    6.  Substituting K(x) into the force equation gives the final result.
    """
    # Define symbolic variables
    x, l, n, E0 = sympy.symbols('x l n E(0)')

    # Build the terms of the equation
    coefficient = (2 * E0 * x) / (n**2 * l**2)
    exponential_term = sympy.exp(x**2 / (n**2 * l**2))

    # The final force equation
    force_equation = coefficient * exponential_term

    # Print the equation and its components
    print("The derived force law F(x) between the polymer ends is:")
    print(f"F(x) = {sympy.pretty(force_equation, use_unicode=False)}")
    print("\nThis can be written as F(x) = (A) * exp(B), where:")
    print(f"  A = {sympy.pretty(coefficient, use_unicode=False)}")
    print(f"  B = {sympy.pretty(exponential_term.exp, use_unicode=False)}")

    # The problem asks to output each number in the final equation.
    # Let's spell out the numeric constants present.
    print("\nThe numeric constants in the final formula are:")
    print("The leading coefficient is: 2")
    print("The exponent of n in the denominator is: 2")
    print("The exponent of l in the denominator is: 2")
    print("The exponent of x in the exponential is: 2")
    print("The exponent of n in the exponential denominator is: 2")
    print("The exponent of l in the exponential denominator is: 2")

solve_polymer_force()