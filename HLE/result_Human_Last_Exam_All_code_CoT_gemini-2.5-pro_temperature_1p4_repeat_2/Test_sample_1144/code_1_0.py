import math

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron
    to its momentum in the first Bohr orbit.
    """

    # --- Physical Constants and Given Values ---

    # Bohr radius (radius of the first Bohr orbit), in meters
    bohr_radius_a0 = 5.29177e-11

    # Given uncertainty in position, in meters (10 pm = 10 * 10^-12 m)
    uncertainty_position_dx = 10e-12

    # --- Calculation ---

    # The ratio of momentum uncertainty to momentum (Δp / p) simplifies to:
    # (ħ / (2 * Δx)) / (ħ / a₀) = a₀ / (2 * Δx)
    # where ħ is the reduced Planck constant.
    # We can directly calculate this simplified ratio.

    ratio = bohr_radius_a0 / (2 * uncertainty_position_dx)

    # --- Output ---
    print("This script calculates the ratio of the uncertainty in momentum (Δp) to the momentum (p) of an electron in the first Bohr orbit.")
    print("The formula for the ratio Δp / p simplifies to a₀ / (2 * Δx).")
    print("\n--- Final Equation ---")
    print(f"Ratio = Bohr Radius / (2 * Uncertainty in Position)")
    print(f"Ratio = {bohr_radius_a0} m / (2 * {uncertainty_position_dx} m)")
    print(f"Final Ratio = {ratio}")

if __name__ == "__main__":
    calculate_momentum_ratio()