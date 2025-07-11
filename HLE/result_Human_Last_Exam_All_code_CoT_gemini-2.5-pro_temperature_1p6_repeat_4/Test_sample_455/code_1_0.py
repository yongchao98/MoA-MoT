import numpy as np
from scipy.special import ai_zeros

def calculate_energy_difference():
    """
    This function calculates the energy difference between the first and second
    energy levels of a quantum particle in a specified potential well.
    
    The calculation follows these steps:
    1.  The complex potential given is simplified. We analyze the potential for r > R,
        V(r) = V0 * (1 - (R/r)**2), and assume an infinite wall at r=R. This is
        justified because the potential for r < R is a large, repulsive barrier.
    2.  The 3D radial Schrödinger equation is set up for this potential. The solutions
        are related to the modified Bessel function of the second kind with imaginary
        order, K_{i*beta}(k*r).
    3.  The boundary condition u(R) = 0 implies that the energy levels are quantized
        by the zeros of this special function. The energy is given by
        E_nl = V0 - (hbar^2 * z_nl^2) / (2*m*R^2), where z_nl is the n-th zero.
    4.  The ground state (E1) corresponds to n=1, l=0. The first excited state (E2)
        is the next highest energy, which is found to be the n=1, l=1 state.
    5.  An asymptotic formula is used to find the required zeros of the Bessel function.
    6.  The energies E1 and E2 are calculated, and their difference is found.
    """
    
    # Constants and Parameters
    V0_eV = 15.0  # Potential parameter in eV
    R = 3.0e-9    # Well radius in meters
    m = 9.11e-31  # Mass of the particle (electron) in kg
    hbar = 1.05457e-34  # Reduced Planck constant in J*s
    e_charge = 1.60218e-19  # Elementary charge in C for eV to J conversion
    
    # Convert V0 to Joules
    V0_J = V0_eV * e_charge
    
    # --- Calculation for the l=0 state (Ground State, E1) ---
    l0 = 0
    # Calculate gamma and beta for l=0
    gamma0 = (2 * m * V0_J * R**2) / hbar**2
    beta0 = np.sqrt(gamma0 - l0 * (l0 + 1) - 0.25)
    
    # --- Calculation for the l=1 state (First Excited State, E2) ---
    l1 = 1
    # Calculate beta for l=1
    beta1 = np.sqrt(gamma0 - l1 * (l1 + 1) - 0.25)

    # Use the first zero of the Airy function Ai(-x) for the lowest radial state (n=1)
    a1 = ai_zeros(1)[0][0] # This gives a_1 ≈ 2.3381
    
    # Asymptotic formula for the first zero (n=1) of K_{i*beta}(z)
    z_1_0 = beta0 - a1 * (2 * beta0)**(1/3)
    z_1_1 = beta1 - a1 * (2 * beta1)**(1/3)

    # Calculate the energy factor hbar^2 / (2*m*R^2)
    energy_factor_J = hbar**2 / (2 * m * R**2)
    
    # Calculate the energies E1 and E2 in Joules
    E1_J = V0_J - energy_factor_J * z_1_0**2
    E2_J = V0_J - energy_factor_J * z_1_1**2

    # Calculate the energy difference in Joules
    delta_E_J = E2_J - E1_J
    
    # Convert results to eV for printing
    E1_eV = E1_J / e_charge
    E2_eV = E2_J / e_charge
    delta_E_eV = delta_E_J / e_charge
    
    # Print the required output
    print(f"The first energy level E1 (n=1, l=0) is: {E1_eV:.4f} eV")
    print(f"The second energy level E2 (n=1, l=1) is: {E2_eV:.4f} eV")
    print(f"The energy difference is \u0394E = E2 - E1 = {E2_eV:.4f} eV - {E1_eV:.4f} eV = {delta_E_eV:.4f} eV")
    
    # Final answer in the required format
    print(f"\n<<<{delta_E_eV:.4f}>>>")

calculate_energy_difference()