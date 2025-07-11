import scipy.constants

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron
    to its momentum in the first Bohr orbit.
    """
    # Known values
    # Position uncertainty given in the problem (10 pm)
    delta_x = 10e-12  # in meters

    # Physical constants
    # Bohr radius (radius of the first Bohr orbit)
    bohr_radius = scipy.constants.physical_constants['Bohr radius'][0] # a_0

    # The ratio Δp/p simplifies to a_0 / (2 * Δx)
    # Δp ≈ ħ / (2 * Δx)
    # p = ħ / a_0 (for the first Bohr orbit, n=1)
    # Ratio = (ħ / (2 * Δx)) / (ħ / a_0) = a_0 / (2 * Δx)
    ratio = bohr_radius / (2 * delta_x)

    # Output the result and the equation with the values used
    print("To find the ratio of the uncertainty of momentum (Δp) to the momentum (p), we use the formula:")
    print("Ratio = a₀ / (2 * Δx)\n")
    print("Where:")
    print(f"  Bohr radius (a₀) = {bohr_radius:.5e} m")
    print(f"  Position uncertainty (Δx) = {delta_x:.5e} m\n")
    print("Calculation:")
    print(f"Ratio = {bohr_radius:.5f} / (2 * {delta_x:.5f})")
    print(f"Ratio = {ratio:.6f}")

if __name__ == "__main__":
    calculate_momentum_ratio()