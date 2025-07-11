import math

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of an electron's momentum to its
    momentum in the first Bohr orbit.
    """
    # Constants and given values
    # Bohr radius (r1) for the first orbit in meters
    r1 = 5.29177210903e-11
    # Uncertainty in position (delta_x) in meters
    # Given as 10 pm = 10 * 10^-12 m
    delta_x = 10e-12

    # The ratio of the uncertainty in momentum (delta_p) to the momentum (p) is
    # derived as:
    # delta_p / p = (hbar / (2 * delta_x)) / (hbar / r1) = r1 / (2 * delta_x)
    ratio = r1 / (2 * delta_x)

    print("This script calculates the ratio of momentum uncertainty to momentum (Δp/p).")
    print("The simplified formula for the ratio is r1 / (2 * Δx).")
    print("-" * 50)
    print(f"Radius of the first Bohr orbit (r1): {r1:.5e} m")
    print(f"Uncertainty in position (Δx): {delta_x:.5e} m")
    print("-" * 50)
    print("Final Calculation:")
    print(f"Ratio = {r1:.5e} / (2 * {delta_x:.5e})")
    print(f"Ratio = {ratio:.6f}")

if __name__ == "__main__":
    calculate_momentum_ratio()