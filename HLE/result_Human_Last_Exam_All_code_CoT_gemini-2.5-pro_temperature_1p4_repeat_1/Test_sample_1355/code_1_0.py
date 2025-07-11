import sympy
from sympy import integrate, pi, sqrt, Symbol

def solve_conductance_ratio():
    """
    Calculates the ratio of the fourth to the first moment of conductance
    for a system in symmetry Class D using Random Matrix Theory.
    """
    # Define the symbolic variable for dimensionless conductance 'g'.
    # In this model, g represents the transmission eigenvalue T.
    g = Symbol('g', real=True)

    # Define the probability distribution P(g) for a single-channel scatterer
    # in symmetry class D. This is the arcsine distribution.
    # P(g) = 1 / (pi * sqrt(g * (1 - g)))
    P_g = 1 / (pi * sqrt(g * (1 - g)))

    # Calculate the average conductance, which is the first statistical moment <g>.
    # <g> = integral from 0 to 1 of g * P(g) dg
    try:
        avg_g = integrate(g * P_g, (g, 0, 1))
    except Exception as e:
        print(f"Error calculating the average conductance: {e}")
        return

    # Calculate the fourth statistical moment of the conductance, <g^4>.
    # <g^4> = integral from 0 to 1 of g^4 * P(g) dg
    try:
        fourth_moment_g = integrate(g**4 * P_g, (g, 0, 1))
    except Exception as e:
        print(f"Error calculating the fourth moment: {e}")
        return

    # Ensure the moments are not zero to avoid division errors.
    if avg_g == 0:
        print("Error: Average conductance is zero, cannot compute the ratio.")
        return

    # Calculate the final ratio.
    ratio = fourth_moment_g / avg_g

    # Print the equation and the final result.
    print("The problem asks for the ratio: <g^4> / <g>")
    print(f"The calculated average value <g> is: {avg_g}")
    print(f"The calculated fourth moment <g^4> is: {fourth_moment_g}")
    print(f"The final ratio (<g^4> / <g>) is: {fourth_moment_g} / {avg_g} = {ratio}")

if __name__ == "__main__":
    solve_conductance_ratio()
