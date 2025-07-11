import math

def solve_energy_loss():
    """
    Calculates the energy loss per centimeter for alpha particles in air
    based on Geiger's rule for range and energy.
    """
    # Given values
    E0 = 8.5  # Initial energy in MeV
    R0 = 8.3  # Total range in cm
    x = 4.0   # Distance from the source in cm

    print(f"Problem: Calculate the energy loss per cm for an alpha particle at x = {x} cm.")
    print(f"Given: Initial Energy (E0) = {E0} MeV, Total Range (R0) = {R0} cm.")
    print("-" * 50)

    # Step 1: Calculate the constant 'a' from Geiger's rule R = a * E^(3/2)
    # a = R0 / E0^(3/2)
    a = R0 / math.pow(E0, 1.5)
    print("Step 1: Find the constant 'a' using Geiger's rule R = a*E^(3/2).")
    print(f"a = R0 / E0^(3/2) = {R0} / {E0}^(1.5) = {a:.4f}")
    print("-" * 50)

    # Step 2: Calculate the remaining range R' and the energy E_x at distance x
    R_prime = R0 - x
    # E_x = (R' / a)^(2/3)
    E_at_x = math.pow(R_prime / a, 2/3)
    print(f"Step 2: Find the particle's energy (E_x) at x = {x} cm.")
    print(f"The remaining range (R') = R0 - x = {R0} - {x} = {R_prime:.1f} cm.")
    print(f"The energy at this point E_x = (R' / a)^(2/3) = ({R_prime:.1f} / {a:.4f})^(2/3) = {E_at_x:.4f} MeV.")
    print("-" * 50)

    # Step 3: Calculate the energy loss per centimeter, dE/dx
    # From R = a*E^(3/2), we can derive dE/dx = (2 / (3*a)) * E^(-1/2)
    dEdx = (2 / (3 * a)) * math.pow(E_at_x, -0.5)
    print("Step 3: Calculate the energy loss per centimeter (dE/dx).")
    print("The formula is dE/dx = (2 / (3*a)) * E_x^(-1/2).")
    
    # Output the final equation with all the numbers, as requested.
    print("\nFinal Equation:")
    print(f"dE/dx = (2 / (3 * {a:.4f})) * {E_at_x:.4f}^(-0.5)")
    print(f"dE/dx = {dEdx:.3f} MeV/cm")
    print("-" * 50)

    print(f"\nFinal Answer: The energy loss per centimetre at a distance of {x} cm is {dEdx:.3f} MeV/cm.")

solve_energy_loss()