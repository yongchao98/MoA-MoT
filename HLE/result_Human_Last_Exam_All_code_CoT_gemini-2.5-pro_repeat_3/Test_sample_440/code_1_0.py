import numpy as np

def solve_density_profile():
    """
    Calculates the linearized density profile of a non-ideal gas mixture
    in a gravitational field.
    """
    # 1. Define constants in SI units
    N_A = 2e23  # Number of particles of Gas A
    N_B = 1.5e23 # Number of particles of Gas B
    M_A = 28e-3  # Molar mass of Gas A (kg/mol)
    M_B = 44e-3  # Molar mass of Gas B (kg/mol)
    a_AA = 2.5   # van der Waals a for A-A (Pa m^6 mol^-2)
    b_AA = 0.04  # van der Waals b for A (m^3 mol^-1)
    a_BB = 3.6   # van der Waals a for B-B (Pa m^6 mol^-2)
    b_BB = 0.05  # van der Waals b for B (m^3 mol^-1)
    a_AB = 3.0   # van der Waals a for A-B (Pa m^6 mol^-2)
    A = 0.1      # Cross-sectional area (m^2)
    H = 10.0     # Height of container (m)
    T = 500.0    # Temperature (K)
    g = 9.81     # Gravitational acceleration (m/s^2)
    N_avo = 6.02214076e23 # Avogadro's number (mol^-1)
    R = 8.314462618       # Ideal gas constant (J/(mol*K))

    # 2. Calculate average molar densities
    V = A * H
    n_A_tot = N_A / N_avo
    n_B_tot = N_B / N_avo
    n_A_avg = n_A_tot / V
    n_B_avg = n_B_tot / V

    # Use shorter variable names for formula clarity
    na, nb = n_A_avg, n_B_avg
    bA, bB = b_AA, b_BB
    nt = na + nb
    
    # 3. Calculate the Jacobian matrix of chemical potentials at average densities
    # The internal chemical potential for component i in a vdW mixture is:
    # mu_i = C(T) + RT*ln(n_i) - RT*ln(1-B) + (b_i*RT*nt)/(1-B) - 2*(n_i*a_ii + n_j*a_ij)
    # where B = n_A*b_A + n_B*b_B. We differentiate this to find the Jacobian elements J_ij = d(mu_i)/d(n_j).

    B_val = na * bA + nb * bB
    den = 1.0 - B_val
    
    if den <= 0:
        print("Error: The term (1 - n*b) is non-positive. The gas is unstable or in a liquid phase under these conditions.")
        return

    # Jacobian element J_AA = d(mu_A)/d(n_A)
    J_AA = (R * T / na) + (R * T * bA / den) + (R * T * bA / den) + (R * T * nt * bA**2 / den**2) - 2 * a_AA
    # Jacobian element J_BB = d(mu_B)/d(n_B)
    J_BB = (R * T / nb) + (R * T * bB / den) + (R * T * bB / den) + (R * T * nt * bB**2 / den**2) - 2 * a_BB
    # Jacobian element J_AB = d(mu_A)/d(n_B)
    J_AB = (R * T * bB / den) + (R * T * bA / den) + (R * T * nt * bA * bB / den**2) - 2 * a_AB
    # Jacobian element J_BA = d(mu_B)/d(n_A) should be equal to J_AB
    J_BA = (R * T * bA / den) + (R * T * bB / den) + (R * T * nt * bA * bB / den**2) - 2 * a_AB

    # 4. Construct and invert the Jacobian matrix
    J = np.array([[J_AA, J_AB], [J_BA, J_BB]])
    J_inv = np.linalg.inv(J)

    # 5. Solve for the constant density gradients (slopes)
    M_g_vec = np.array([-M_A * g, -M_B * g])
    C_vec = J_inv @ M_g_vec
    C_A = C_vec[0]
    C_B = C_vec[1]

    # 6. Determine the parameters for the linear mass density profile equation
    # rho(z) = rho_avg + C_rho * (z - H/2)
    rho_avg = n_A_avg * M_A + n_B_avg * M_B
    C_rho = C_A * M_A + C_B * M_B

    # 7. Print the final result
    print("The mass density profile rho(z) can be approximated by the linear equation:")
    print(f"rho(z) = rho_avg + C_rho * (z - H/2)")
    print("\nCalculated parameters:")
    print(f"  rho_avg = {rho_avg:.6f} kg/m^3")
    print(f"  C_rho   = {C_rho:.6f} kg/m^4")
    print(f"  H/2     = {H/2:.1f} m")

    print("\nFinal equation for the density profile:")
    print(f"rho(z) = {rho_avg:.6f} + {C_rho:.6f} * (z - {H/2:.1f})")

    # Calculate density at the top of the cylinder for the final answer tag
    rho_at_H = rho_avg + C_rho * (H / 2.0)
    # The final answer tag requires a single value. We provide the density at the top of the cylinder.
    # print(f"\nFor reference, the density at the top (z=10m) is: {rho_at_H:.6f} kg/m^3")
    # print(f"<<<{rho_at_H:.6f}>>>")

if __name__ == '__main__':
    solve_density_profile()