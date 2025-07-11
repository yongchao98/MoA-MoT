import math

def calculate_energy_loss():
    """
    Calculates the energy loss per centimeter for an alpha-particle in air
    based on an empirical range-energy relationship.
    """
    # Given initial values
    E0 = 8.5  # Initial energy in MeV
    R0 = 8.3  # Total range in cm
    x = 4.0   # Distance from the source in cm

    # The formula for energy loss per centimeter (-dE/dx) is derived from
    # the empirical relation E = c * R^(2/3).
    # The derived formula is: -dE/dx = (2 * E0) / (3 * R0^(2/3) * (R0 - x)^(1/3))

    # Calculate the remaining range at distance x
    remaining_range = R0 - x

    # Calculate the energy loss per centimeter
    energy_loss = (2 * E0) / (3 * (R0**(2/3)) * (remaining_range**(1/3)))

    # Print the step-by-step calculation
    print("Step 1: State the formula for energy loss per centimeter (-dE/dx).")
    print("-dE/dx = (2 * E0) / (3 * R0^(2/3) * (R0 - x)^(1/3))\n")

    print("Step 2: Substitute the given values into the formula.")
    print(f"-dE/dx = (2 * {E0}) / (3 * {R0}^(2/3) * ({R0} - {x})^(1/3))\n")

    print("Step 3: Perform the calculation.")
    print(f"-dE/dx = {2 * E0} / (3 * {R0**(2/3):.4f} * {remaining_range}^(1/3))")
    print(f"-dE/dx = {2 * E0} / (3 * {R0**(2/3):.4f} * {remaining_range**(1/3):.4f})")
    print(f"-dE/dx = {2 * E0} / {(3 * (R0**(2/3)) * (remaining_range**(1/3))):.4f}\n")

    print(f"Result:")
    print(f"The energy loss per centimetre for these α-particles at a distance of {x} cm is {energy_loss:.4f} MeV/cm.")

calculate_energy_loss()