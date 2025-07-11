import math

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron
    to its momentum in the first Bohr orbit.
    """
    # Constants
    # Bohr radius (r1) in meters for the first orbit (n=1).
    r1 = 5.29177e-11  # meters

    # Given uncertainty in position in picometers.
    delta_x_pm = 10.0

    # Convert the uncertainty in position from picometers to meters.
    delta_x_m = delta_x_pm * 1e-12  # meters

    # The ratio is derived as r1 / (2 * delta_x).
    # delta_p = h_bar / (2 * delta_x) from Heisenberg's principle.
    # p = h_bar / r1 from Bohr's model for the first orbit.
    # The ratio delta_p / p simplifies to r1 / (2 * delta_x).
    ratio = r1 / (2 * delta_x_m)

    # Print the final equation with the numbers plugged in.
    print("The ratio is calculated using the formula: r1 / (2 * delta_x)")
    print("Where:")
    print(f"  r1 (Bohr Radius) = {r1} m")
    print(f"  delta_x (Position Uncertainty) = {delta_x_m} m")
    print("\nCalculation:")
    print(f"Ratio = {r1} / (2 * {delta_x_m})")
    print(f"\nThe calculated ratio is: {ratio}")

if __name__ == "__main__":
    calculate_momentum_ratio()