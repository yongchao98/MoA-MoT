import math

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron
    to its momentum in the first Bohr orbit.
    """
    # Bohr radius (a0) in meters
    a0 = 5.29177210903e-11

    # Given uncertainty in position (delta_x) in picometers, converted to meters
    delta_x_pm = 10
    delta_x_m = delta_x_pm * 1e-12

    # The ratio is derived from Heisenberg's principle and the Bohr model:
    # Uncertainty in momentum, delta_p = hbar / (2 * delta_x)
    # Momentum in first Bohr orbit, p = hbar / a0
    # Ratio = delta_p / p = (hbar / (2 * delta_x)) / (hbar / a0) = a0 / (2 * delta_x)
    ratio = a0 / (2 * delta_x_m)

    # Print the equation with the values plugged in
    print("This script calculates the ratio of momentum uncertainty to momentum for an electron in the first Bohr orbit.")
    print("Formula: Ratio = Bohr Radius / (2 * Position Uncertainty)")
    print(f"Ratio = {a0} m / (2 * {delta_x_m} m)")
    print(f"The calculated ratio is: {ratio}")

if __name__ == "__main__":
    calculate_momentum_ratio()