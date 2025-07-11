import math

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron
    to its momentum in the first Bohr orbit.
    """
    # Physical constants
    h_bar = 1.054571817e-34  # Reduced Planck constant (J*s)
    bohr_radius = 5.29177210903e-11  # Radius of the first Bohr orbit (m)

    # Given value
    delta_x = 10e-12  # Uncertainty in position (10 pm in meters)

    # 1. Calculate the uncertainty in momentum (Δp)
    delta_p = h_bar / (2 * delta_x)

    # 2. Calculate the momentum of the electron in the first Bohr orbit (p1)
    p1 = h_bar / bohr_radius

    # 3. Calculate the ratio
    ratio = delta_p / p1

    # Print the results in the requested format
    print(f"Given uncertainty in position (Δx) = {delta_x:.2e} m")
    print(f"Bohr radius (r₁) = {bohr_radius:.2e} m")
    print("-" * 30)
    print(f"Uncertainty in momentum (Δp) = {delta_p:.4e} kg*m/s")
    print(f"Momentum in first Bohr orbit (p₁) = {p1:.4e} kg*m/s")
    print("-" * 30)
    print("Ratio (Δp / p₁):")
    print(f"{delta_p:.4e} / {p1:.4e} = {ratio:.4f}")

if __name__ == "__main__":
    calculate_momentum_ratio()