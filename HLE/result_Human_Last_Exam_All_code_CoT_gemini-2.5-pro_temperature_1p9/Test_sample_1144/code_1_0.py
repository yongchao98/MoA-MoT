import math

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron to its
    momentum in the first Bohr orbit.
    """
    # 1. Define the known physical constants and given values.
    # Radius of the first Bohr orbit (r_1) in meters.
    r_1 = 5.29177e-11
    # Given uncertainty in position (delta_x) is 10 pm, which is 10e-12 meters.
    delta_x = 10e-12

    # 2. Calculate the ratio using the simplified formula derived from the principles:
    # Ratio = (delta_p / p_1) = r_1 / (2 * delta_x)
    ratio = r_1 / (2 * delta_x)

    # 3. Print the final result and the equation with the numbers plugged in.
    print("The problem is to find the ratio of momentum uncertainty to momentum (Δp/p).")
    print("This simplifies to the formula: Ratio = r₁ / (2 * Δx)")
    print("\nGiven values:")
    print(f"Bohr Radius (r₁): {r_1} m")
    print(f"Position Uncertainty (Δx): {delta_x} m")
    print("\nCalculation:")
    print(f"Ratio = {r_1} / (2 * {delta_x})")
    print(f"Final Ratio = {ratio}")

if __name__ == "__main__":
    calculate_momentum_ratio()