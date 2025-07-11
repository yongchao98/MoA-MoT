import sympy
from sympy import integrate, Symbol, pi, sqrt, S

def solve_conductance_ratio():
    """
    Calculates the ratio of the fourth moment to the first moment of conductance
    for a disordered Majorana wire at the critical point.
    """
    # Define the symbol for the dimensionless conductance g
    g = Symbol('g')

    # Define the arcsine probability distribution function for g
    # P(g) = 1 / (pi * sqrt(g * (1 - g)))
    P_g = 1 / (pi * sqrt(g * (1 - g)))

    print("This script calculates the ratio <g^4> / <g> for the dimensionless conductance `g` of a Majorana wire at the critical point.")
    print("The conductance follows the arcsine distribution: P(g) = 1 / (pi * sqrt(g*(1-g))).")
    print("-" * 50)

    # --- Step 1: Calculate the average conductance <g> ---
    print("Step 1: Calculating the average conductance <g>.")
    # The integrand is g * P(g)
    integrand_g1 = g * P_g
    # Perform the definite integral from g = 0 to g = 1
    average_g = integrate(integrand_g1, (g, 0, 1))
    print(f"The average value is <g> = Integral(g * P(g) dg) from 0 to 1.")
    print(f"Result: <g> = {average_g}")
    print("-" * 50)

    # --- Step 2: Calculate the fourth statistical moment <g^4> ---
    print("Step 2: Calculating the fourth statistical moment <g^4>.")
    # The integrand is g^4 * P(g)
    integrand_g4 = g**4 * P_g
    # Perform the definite integral from g = 0 to g = 1
    fourth_moment_g = integrate(integrand_g4, (g, 0, 1))
    print(f"The fourth moment is <g^4> = Integral(g^4 * P(g) dg) from 0 to 1.")
    print(f"Result: <g^4> = {fourth_moment_g}")
    print("-" * 50)

    # --- Step 3: Compute the ratio <g^4> / <g> ---
    print("Step 3: Computing the final ratio.")
    if average_g.is_zero:
        ratio = "undefined (average conductance is zero)"
    else:
        # Calculate the ratio
        ratio = fourth_moment_g / average_g
    
    print("The final equation is: <g^4> / <g>")
    # Using S() to ensure sympy objects are used for pretty printing
    final_equation = f"{S(fourth_moment_g)} / {S(average_g)} = {S(ratio)}"
    print(final_equation)

if __name__ == '__main__':
    solve_conductance_ratio()
