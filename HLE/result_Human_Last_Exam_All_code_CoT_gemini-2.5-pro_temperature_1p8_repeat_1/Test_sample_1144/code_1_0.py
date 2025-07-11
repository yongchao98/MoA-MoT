import math

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron
    to its momentum in the first Bohr orbit.
    """
    # Step 1: Define physical constants and given values.
    # Reduced Planck constant in J·s
    hbar = 1.054571817e-34
    # Bohr radius (radius of the first Bohr orbit) in meters
    a0 = 5.29177210903e-11
    # Given uncertainty in position in meters (10 pm = 10 * 10^-12 m)
    delta_x = 10e-12

    # Step 2: Calculate the uncertainty in momentum (delta_p)
    # using Heisenberg's Uncertainty Principle (Δp * Δx ≥ ħ/2).
    # We calculate the minimum uncertainty.
    delta_p = hbar / (2 * delta_x)

    # Step 3: Calculate the momentum (p) of the electron in the first Bohr orbit.
    # From Bohr's model, angular momentum L = p * r = n * ħ.
    # For the first orbit (n=1), p = ħ / a₀.
    p_bohr = hbar / a0

    # Step 4: Calculate the final ratio of Δp to p.
    ratio = delta_p / p_bohr

    # Print the final result, showing each number in the final equation as requested.
    print(f"The uncertainty in momentum (Δp) is: {delta_p:.5g} kg·m/s")
    print(f"The momentum in the first Bohr orbit (p) is: {p_bohr:.5g} kg·m/s")
    print("\nThe ratio (Δp / p) is calculated as follows:")
    print(f"{delta_p} / {p_bohr} = {ratio}")

if __name__ == "__main__":
    calculate_momentum_ratio()