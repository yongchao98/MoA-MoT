import scipy.constants

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron
    to its momentum in the first Bohr orbit.
    """
    # 1. Define constants and given values
    # Reduced Planck constant in J·s
    hbar = scipy.constants.hbar
    # Bohr radius (radius of the first Bohr orbit) in meters
    a0 = scipy.constants.physical_constants['Bohr radius'][0]
    # Given uncertainty in position in meters (10 pm = 10e-12 m)
    delta_x = 10e-12

    # 2. Calculate the uncertainty in momentum (delta_p) using Heisenberg's Uncertainty Principle
    # Δp = ħ / (2 * Δx)
    delta_p = hbar / (2 * delta_x)

    # 3. Calculate the momentum (p) in the first Bohr orbit
    # p = ħ / a₀
    p = hbar / a0

    # 4. Calculate the final ratio
    ratio = delta_p / p

    # 5. Print the results in a clear, step-by-step format
    print("Calculating the ratio (Δp / p):")
    print("-" * 40)
    print(f"Uncertainty in momentum (Δp) = ħ / (2 * Δx)")
    print(f"Δp = {hbar:.5e} / (2 * {delta_x:.5e}) = {delta_p:.5e} kg·m/s")
    print("-" * 40)
    print(f"Momentum in first Bohr orbit (p) = ħ / a₀")
    print(f"p = {hbar:.5e} / {a0:.5e} = {p:.5e} kg·m/s")
    print("-" * 40)
    print("Final Ratio = Δp / p")
    # The final equation with all the numbers
    print(f"Ratio = {delta_p:.5e} / {p:.5e} = {ratio:.4f}")

if __name__ == "__main__":
    calculate_momentum_ratio()