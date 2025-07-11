import numpy as np

def calculate_lindhard_function_static_limit():
    """
    Calculates the Lindhard polarization function for a 3D homogeneous electron gas
    at T=0 in the static, long-wavelength limit (q=0, omega=0).

    This value is equal to the density of states at the Fermi energy, g(epsilon_F).
    The calculation is performed for Sodium (Na) as a representative metal.
    """

    # Physical constants in SI units
    m_e = 9.1093837e-31  # Electron mass (kg)
    hbar = 1.054571817e-34 # Reduced Planck constant (J*s)

    # Electron number density for Sodium (Na) in m^-3
    n = 2.54e28

    print("The Lindhard function at q=0, omega=0 is equal to the density of states at the Fermi energy g(eps_F).")
    print("The electron density 'n' for Sodium is used for a numerical example.")
    print(f"n = {n:.3e} m^-3\n")

    # Step 1: Calculate the Fermi wavevector k_F
    # Formula: k_F = (3 * pi^2 * n)^(1/3)
    k_f = (3 * np.pi**2 * n)**(1/3)
    print(f"Step 1: Calculate the Fermi wavevector, k_F.")
    print(f"k_F = (3 * pi^2 * n)^(1/3) = (3 * {np.pi**2:.4f} * {n:.3e})^(1/3) = {k_f:.4e} m^-1\n")

    # Step 2: Calculate the Lindhard function Pi_0(0,0) = g(eps_F)
    # Formula: g(eps_F) = (m_e * k_F) / (pi^2 * hbar^2)
    lindhard_value = (m_e * k_f) / (np.pi**2 * hbar**2)

    print("Step 2: Calculate the numerical value using the formula Pi_0(0,0) = g(eps_F) = (m_e * k_F) / (pi^2 * hbar^2).")
    print("Substituting the numerical values into the equation:")
    print(f"Pi_0(0,0) = ({m_e:.4e} kg * {k_f:.4e} m^-1) / ({np.pi**2:.4f} * ({hbar:.4e} J*s)^2)")
    print(f"Pi_0(0,0) = {m_e * k_f:.4e} / {np.pi**2 * hbar**2:.4e}")
    print(f"The final calculated value is: {lindhard_value:.4e} J^-1 m^-3")

    # The result in more common condensed matter units (eV^-1 Angstrom^-3)
    eV = 1.60218e-19 # Joules
    angstrom = 1e-10 # meters
    lindhard_value_ev_a = lindhard_value * eV * (angstrom**3)
    # print(f"(For reference, this is ~{lindhard_value_ev_a:.3f} eV^-1 Angstrom^-3)")
    
    # Return final answer for the '<<<' block
    return lindhard_value

if __name__ == '__main__':
    final_answer = calculate_lindhard_function_static_limit()
    # print(f"\n<<<final_answer>>>")

calculate_lindhard_function_static_limit()