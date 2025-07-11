import numpy as np
from scipy.integrate import quad

def solve_conductance_ratio():
    """
    Calculates the ratio of the fourth to the first statistical moment of
    conductance for disordered Majorana wires at the critical point.

    The problem reduces to calculating moments of the conductance g, whose
    probability distribution P(g) is known from random matrix theory for
    symmetry class D.
    P(g) = 1 / (2 * sqrt(g)) for g in (0, 1].

    The n-th moment <g^n> is the integral of g^n * P(g) from 0 to 1.
    """

    # The integrand for the average value <g> = <g^1>
    # The expression is g * (1 / (2 * sqrt(g))) = sqrt(g) / 2
    def integrand_g1(g):
        return g**0.5 / 2

    # The integrand for the fourth moment <g^4>
    # The expression is g^4 * (1 / (2 * sqrt(g))) = g^3.5 / 2
    def integrand_g4(g):
        return g**3.5 / 2

    # --- Calculation ---
    # The integration range is from 0 to 1.
    # quad() returns the result and an estimated error. We only need the result.
    
    # Calculate the average value <g>
    avg_g, _ = quad(integrand_g1, 0, 1)

    # Calculate the fourth moment <g^4>
    moment_g4, _ = quad(integrand_g4, 0, 1)

    # Calculate the final ratio
    if avg_g == 0:
        ratio = float('inf')
    else:
        ratio = moment_g4 / avg_g

    # --- Output ---
    print("Calculating the ratio between the fourth statistical moment and the average value of conductance.")
    print("-" * 80)
    print(f"The average conductance <g> is the integral of g * P(g):")
    print(f"<g> = {avg_g:.8f}")
    print("-" * 80)
    print(f"The fourth moment of conductance <g^4> is the integral of g^4 * P(g):")
    print(f"<g^4> = {moment_g4:.8f}")
    print("-" * 80)
    print("The final ratio is <g^4> / <g>:")
    print(f"{moment_g4:.8f} / {avg_g:.8f} = {ratio:.8f}")
    print("-" * 80)
    print(f"The exact analytical result is (1/9) / (1/3) = 1/3.")

solve_conductance_ratio()