import numpy as np

def solve_density_profile():
    """
    Calculates and prints the density profile of a two-gas mixture in a
    vertical cylinder under gravity.
    """
    # Given Parameters (in SI units)
    N_A = 2e23       # Number of particles of Gas A
    N_B = 1.5e23     # Number of particles of Gas B
    M_A = 0.028      # Molar mass of Gas A (kg/mol)
    M_B = 0.044      # Molar mass of Gas B (kg/mol)
    A = 0.1          # Cross-sectional area of the cylinder (m^2)
    H = 10.0         # Height of the container (m)
    T = 500.0        # Temperature (K)
    g = 9.81         # Gravitational acceleration (m/s^2)

    # Physical Constants
    R = 8.31446      # Ideal gas constant (J/(mol*K))
    N_AVOGADRO = 6.02214e23 # Avogadro's number (mol^-1)

    # The problem states the system obeys the Boltzmann distribution in the
    # gravitational field. This means the density profile is determined by the
    # barometric formula, which depends on gravitational potential energy.
    # The van der Waals parameters affect the equation of state (pressure) but
    # are not required for calculating the density profile in this context.

    # Step 1: Calculate the exponential decay constants k_i = M_i * g / (R * T)
    k_A = (M_A * g) / (R * T)
    k_B = (M_B * g) / (R * T)

    # Step 2: Determine the number density at z=0 for each gas.
    # The total number of particles N_i is the integral of number density over the volume:
    # N_i = integral from 0 to H of (n_i(0) * exp(-k_i * z) * A) dz
    # Solving this integral for n_i(0) gives:
    # n_i(0) = (N_i * k_i) / (A * (1 - exp(-k_i * H)))

    n_A_0 = (N_A * k_A) / (A * (1 - np.exp(-k_A * H)))
    n_B_0 = (N_B * k_B) / (A * (1 - np.exp(-k_B * H)))

    # Step 3: Convert number density at z=0 to mass density at z=0 (rho_i(0)).
    # rho_i(0) = n_i(0) * mass_of_one_particle = n_i(0) * (M_i / N_AVOGADRO)
    rho_A_0 = n_A_0 * (M_A / N_AVOGADRO)
    rho_B_0 = n_B_0 * (M_B / N_AVOGADRO)

    # Step 4: Print the final equation for the total density profile rho(z).
    # The equation is: rho(z) = rho_A(z) + rho_B(z)
    # rho(z) = rho_A_0 * exp(-k_A * z) + rho_B_0 * exp(-k_B * z)
    print("The density profile rho(z) of the gas mixture is given by the equation:")
    print(f"rho(z) = {rho_A_0:.4f} * exp(-{k_A:.6f} * z) + {rho_B_0:.4f} * exp(-{k_B:.6f} * z)")
    print("\nwhere:")
    print(f"  - The first term represents Gas A.")
    print(f"    - Density at z=0: {rho_A_0:.4f} kg/m^3")
    print(f"    - Decay constant: {k_A:.6f} m^-1")
    print(f"  - The second term represents Gas B.")
    print(f"    - Density at z=0: {rho_B_0:.4f} kg/m^3")
    print(f"    - Decay constant: {k_B:.6f} m^-1")
    print("  - z is the height in meters (m).")
    print("  - rho(z) is the total mass density in kg/m^3.")

    # Format the final answer string for the <<<...>>> block
    global final_answer
    final_answer = f"rho(z) = {rho_A_0:.4f} * exp(-{k_A:.6f} * z) + {rho_B_0:.4f} * exp(-{k_B:.6f} * z)"


# Run the function to get the answer
solve_density_profile()
# The final answer is wrapped in <<<...>>> as requested.
# Using a global variable to pass the result out of the function scope.
print(f"\n<<<{final_answer}>>>")