import math

def solve_energy_loss():
    """
    Calculates the energy loss per centimeter for α-particles in air
    based on Geiger's rule.
    """
    # Given initial values
    E0 = 8.5  # Initial energy in MeV
    R0 = 8.3  # Total range in cm
    x = 4.0   # Distance from the source in cm

    print("This script calculates the energy loss per centimeter for α-particles in standard air.")
    print("We will use the empirical Geiger's rule for the range-energy relationship: R = a * E^(3/2)\n")

    # Step 1: Determine the constant 'a' from the initial conditions.
    print(f"Step 1: Calculate the proportionality constant 'a' using the initial conditions.")
    print(f"Given initial energy E0 = {E0} MeV and total range R0 = {R0} cm.")
    # From R = a * E^(3/2), we get a = R / E^(3/2)
    a = R0 / (E0 ** 1.5)
    print(f"The constant 'a' is calculated as: a = {R0} / {E0}^(3/2) = {a:.4f} cm/MeV^(3/2)\n")

    # Step 2: Find the energy of the α-particle after it has traveled a distance x.
    print(f"Step 2: Find the energy of the α-particle at a distance of x = {x} cm.")
    # The remaining range R_x at distance x is R0 - x.
    R_x = R0 - x
    print(f"The remaining range at {x} cm is R_x = {R0} cm - {x} cm = {R_x} cm.")

    # The energy E_x at this point is found by rearranging Geiger's rule: E = (R/a)^(2/3).
    E_x = (R_x / a) ** (2.0 / 3.0)
    print(f"The corresponding energy at this point is E_x = ({R_x:.1f} / {a:.4f})^(2/3) = {E_x:.3f} MeV\n")

    # Step 3: Calculate the energy loss per centimeter, |dE/dx|.
    print("Step 3: Calculate the energy loss per centimeter (|dE/dx|).")
    print("From R = a * E^(3/2), we find dR/dE = a * (3/2) * E^(1/2).")
    print("The energy loss rate is |dE/dx| = 1 / (dR/dE) = 1 / (a * 1.5 * E^(1/2)).")

    # Calculate dE/dx at the energy E_x.
    dE_dx = 1.0 / (a * 1.5 * (E_x ** 0.5))
    print(f"\nUsing the energy at {x} cm (E_x = {E_x:.3f} MeV), the final calculation is:")
    print(f"Energy Loss = 1 / ({a:.4f} * 1.5 * {E_x:.3f}^(0.5))")
    print(f"The energy loss per centimetre at {x} cm is {dE_dx:.3f} MeV/cm.")

if __name__ == "__main__":
    solve_energy_loss()