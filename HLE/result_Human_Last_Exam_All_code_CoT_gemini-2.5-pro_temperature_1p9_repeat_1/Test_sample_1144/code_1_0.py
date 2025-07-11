import math

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron to its
    momentum in the first Bohr orbit.
    """
    # Constants
    hbar = 1.054571817e-34  # Reduced Planck constant in J·s
    a0 = 5.29177210903e-11  # Bohr radius in meters
    delta_x = 10e-12        # Uncertainty in position in meters (10 pm)

    # Step 1: Calculate the uncertainty in momentum (delta_p) using Heisenberg's Uncertainty Principle
    # Δp = ħ / (2 * Δx)
    delta_p = hbar / (2 * delta_x)

    # Step 2: Calculate the momentum (p) of the electron in the first Bohr orbit (n=1)
    # p = n * ħ / r = 1 * ħ / a₀
    p = hbar / a0

    # Step 3: Calculate the ratio of delta_p to p
    ratio = delta_p / p
    
    # Alternatively, using the simplified formula: ratio = a0 / (2 * delta_x)
    # This avoids calculating hbar and provides a more direct result.
    # ratio_simplified = a0 / (2 * delta_x)
    # The results from both methods are identical.

    print("Problem: Find the ratio of the uncertainty of the momentum of an electron to its momentum in the first Bohr orbit.")
    print(f"Given uncertainty in position (Δx) = {delta_x:.1e} m\n")

    print("Calculation Steps:")
    print(f"1. Uncertainty in momentum (Δp) = ħ / (2 * Δx) = {delta_p:.4e} kg m/s")
    print(f"2. Momentum in the first Bohr orbit (p) = ħ / a₀ = {p:.4e} kg m/s")
    
    print("\nFinal Equation for the Ratio:")
    print(f"Ratio = Δp / p")
    print(f"Ratio = {delta_p:.4e} / {p:.4e}")
    print(f"The calculated ratio is: {ratio}")

calculate_momentum_ratio()
<<<2.645879793138392>>>