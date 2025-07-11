import math

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of an electron's momentum to its momentum
    in the first Bohr orbit.
    """
    # Constants
    # Bohr radius (r1) in meters
    bohr_radius = 5.29177e-11
    # Reduced Planck constant (h-bar) in J·s
    hbar = 1.05457e-34
    # Given uncertainty in position (delta_x) in picometers
    delta_x_pm = 10
    # Convert delta_x from picometers to meters
    delta_x_m = delta_x_pm * 1e-12

    # Step 1: Calculate the uncertainty in momentum (delta_p) using Heisenberg's Uncertainty Principle
    # delta_p * delta_x >= hbar / 2
    # We use the minimum uncertainty: delta_p = hbar / (2 * delta_x)
    delta_p = hbar / (2 * delta_x_m)

    # Step 2: Calculate the momentum (p) of the electron in the first Bohr orbit
    # From Bohr's model, Angular Momentum L = p * r = n * hbar. For n=1, p = hbar / r1.
    p = hbar / bohr_radius

    # Step 3: Calculate the ratio of the uncertainty in momentum to the momentum
    ratio = delta_p / p
    
    # As a simplified check: ratio = (hbar / (2 * delta_x)) / (hbar / bohr_radius) = bohr_radius / (2 * delta_x)
    # simplified_ratio = bohr_radius / (2 * delta_x_m) -> this should be equal to the 'ratio'

    print("Problem: Find the ratio of the uncertainty of the momentum of an electron (Δp) to its momentum (p) in the first Bohr orbit.")
    print(f"Given position uncertainty (Δx) = {delta_x_pm} pm = {delta_x_m:.2e} m")
    print(f"Bohr radius (r₁) = {bohr_radius:.5e} m\n")

    print("Calculations:")
    print(f"1. Uncertainty in momentum (Δp) = ħ / (2 * Δx) = {delta_p:.4e} kg·m/s")
    print(f"2. Momentum in the first Bohr orbit (p) = ħ / r₁ = {p:.4e} kg·m/s")
    
    print("\nFinal Ratio:")
    print(f"The ratio (Δp / p) is {delta_p:.4e} / {p:.4e} = {ratio:.4f}")

calculate_momentum_ratio()
<<<2.6459>>>