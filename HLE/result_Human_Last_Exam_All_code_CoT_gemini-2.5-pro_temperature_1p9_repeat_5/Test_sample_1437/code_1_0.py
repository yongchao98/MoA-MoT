import numpy as np

def calculate_lindhard_zero_limit(electron_density):
    """
    Calculates the Lindhard polarization function at k=0, w=0 for a 3D electron gas.

    This value is equal to the negative of the density of states at the Fermi level, -D(eps_F).
    The calculation is performed for a given electron density in SI units.

    Args:
        electron_density (float): The electron density 'n' in units of m^-3.

    Returns:
        float: The value of the Lindhard function Pi(0,0) in units of J^-1 m^-3.
    """
    # Physical constants in SI units
    m_e = 9.1093837015e-31  # Electron mass in kg
    hbar = 1.054571817e-34    # Reduced Planck constant in J*s
    pi = np.pi

    # Step 1: Calculate the Fermi wavevector (k_F) from the electron density (n)
    # k_F = (3 * pi^2 * n)^(1/3)
    k_F = (3 * pi**2 * electron_density)**(1/3)

    # Step 2: Calculate the density of states at the Fermi level (D(eps_F))
    # D(eps_F) = (m_e * k_F) / (pi^2 * hbar^2)
    dos_fermi_level = (m_e * k_F) / (pi**2 * hbar**2)

    # Step 3: The Lindhard function at k=0, w=0 is -D(eps_F)
    lindhard_value = -dos_fermi_level

    print(f"For an electron density n = {electron_density:.3e} m^-3:")
    print(f"The Fermi wavevector k_F is {k_F:.3e} m^-1.")
    print(f"The density of states at the Fermi level D(eps_F) is {dos_fermi_level:.3e} J^-1 m^-3.")
    print(f"The Lindhard function value Pi(k=0, w=0) is -D(eps_F).")
    print(f"Calculated Numerical Value: {lindhard_value:.3e} J^-1 m^-3")
    
    return lindhard_value

if __name__ == '__main__':
    # Electron density of Sodium (Na), a good approximation for a free electron gas.
    # The Wigner-Seitz radius for Na is r_s = 3.93 a_0, where a_0 is the Bohr radius.
    # n = 3 / (4 * pi * (r_s * a_0)^3)
    a_0 = 5.29177e-11 # Bohr radius in m
    r_s_Na = 3.93
    n_Na = 3 / (4 * pi * (r_s_Na * a_0)**3)
    
    # Calculate and print the result for Sodium
    final_value = calculate_lindhard_zero_limit(n_Na)