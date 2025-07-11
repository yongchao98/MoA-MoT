import scipy.constants

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron
    to its momentum in the first Bohr orbit.
    """
    # Given uncertainty in position in picometers
    delta_x_pm = 10.0

    # Convert picometers to meters
    delta_x_m = delta_x_pm * 1e-12

    # Physical constants from scipy
    hbar = scipy.constants.hbar  # Reduced Planck constant in J·s
    a0 = scipy.constants.physical_constants['Bohr radius'][0]  # Bohr radius in m

    # 1. Calculate the uncertainty in momentum (Δp) using Heisenberg's Uncertainty Principle
    # Δp = ħ / (2 * Δx)
    delta_p = hbar / (2 * delta_x_m)

    # 2. Calculate the momentum of the electron in the first Bohr orbit (p1)
    # For the first orbit (n=1), the radius r1 is the Bohr radius a0.
    # p1 = ħ / r1 = ħ / a0
    p1 = hbar / a0

    # 3. Calculate the ratio of Δp to p1
    ratio = delta_p / p1

    print(f"Given position uncertainty (Δx): {delta_x_pm} pm = {delta_x_m:.2e} m")
    print(f"Bohr radius (a0): {a0:.5e} m\n")
    print(f"Uncertainty in momentum (Δp) = ħ / (2 * Δx) = {delta_p:.5e} kg·m/s")
    print(f"Momentum in first Bohr orbit (p1) = ħ / a0 = {p1:.5e} kg·m/s\n")
    print("Final Ratio Calculation:")
    print(f"Ratio = Δp / p1 = {delta_p:.5e} / {p1:.5e} = {ratio}")

    # The final answer in the required format
    # The prompt requires this specific format at the very end.
    # Note: Using f-string formatting to create the final answer string.
    final_answer_string = f"<<<{ratio}>>>"
    print(final_answer_string)


if __name__ == '__main__':
    calculate_momentum_ratio()