import math

def solve_energy_difference():
    """
    Calculates the energy difference between the first and second energy levels
    for a particle in a 3D infinite spherical potential well.
    """
    # Define physical constants in SI units
    H_BAR = 1.054571817e-34  # Reduced Planck constant in J·s
    M = 9.1093837e-31        # Mass of the particle (electron) in kg
    R = 3e-9                 # Radius of the well in m
    E_CHARGE = 1.602176634e-19 # Elementary charge in C for J-to-eV conversion

    # Zeros of the spherical Bessel functions (z_nl)
    # z_10 is the 1st zero of j_0(x), which is π
    z_10 = math.pi
    # z_11 is the 1st zero of j_1(x)
    z_11 = 4.4934094579

    # The energy levels are given by E_nl = (ħ² * z_nl²) / (2 * m * R²)
    # We can calculate the common energy factor first.
    energy_factor_J = (H_BAR**2) / (2 * M * R**2)

    # Calculate the first energy level (ground state, n=1, l=0) in Joules
    e1_J = energy_factor_J * (z_10**2)

    # Calculate the second energy level (first excited state, n=1, l=1) in Joules
    e2_J = energy_factor_J * (z_11**2)

    # Convert the energy values from Joules to electron volts (eV)
    e1_eV = e1_J / E_CHARGE
    e2_eV = e2_J / E_CHARGE

    # Calculate the energy difference in eV
    delta_e_eV = e2_eV - e1_eV

    # Print the results, including the numbers used in the final equation
    print("Using the infinite spherical well approximation:")
    print(f"The first energy level (ground state) E1 = {e1_eV:.4f} eV")
    print(f"The second energy level (first excited state) E2 = {e2_eV:.4f} eV")
    print("\nThe energy difference ΔE = E2 - E1 is:")
    print(f"{e2_eV:.4f} eV - {e1_eV:.4f} eV = {delta_e_eV:.4f} eV")
    
    return delta_e_eV

# Execute the function and capture the final answer
final_answer = solve_energy_difference()
# print(f"\n<<<{final_answer:.4f}>>>") # Final answer block