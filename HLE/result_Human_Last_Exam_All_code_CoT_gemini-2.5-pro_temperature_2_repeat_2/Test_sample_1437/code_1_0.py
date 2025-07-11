import numpy as np

def calculate_lindhard_function_zero_limit():
    """
    Calculates the numerical value of the Lindhard polarization function
    at zero frequency and zero momentum transfer for a 3D homogeneous
    electron gas, assuming r_s = 1 in Hartree atomic units.
    """
    
    # Step 1: Define constants in Hartree Atomic Units.
    # In this system, h_bar = 1, m_e = 1, e = 1, and the Bohr radius a_0 = 1.
    h_bar = 1.0
    m_e = 1.0
    pi = np.pi

    # Step 2: Set the electron density parameter r_s.
    # A natural choice for a calculation without a specified material is r_s = 1.
    r_s = 1.0
    print(f"Using a jellium model with density parameter r_s = {r_s:.1f} a.u.")

    # Step 3: Calculate electron density 'n' from r_s in atomic units.
    # The volume per electron is (4/3)*pi*(r_s*a_0)^3. With a_0 = 1, n = 1/volume.
    n = 3.0 / (4.0 * pi * r_s**3)

    # Step 4: Calculate the Fermi wavevector k_F from the density n.
    # The relation is n = k_F^3 / (3 * pi^2).
    k_F = (3 * pi**2 * n)**(1.0/3.0)

    # Step 5: The Lindhard function at q=0, omega=0 is -D(epsilon_F).
    # D(epsilon_F) is the density of states at the Fermi level, including spin.
    # The formula is D(epsilon_F) = m_e * k_F / (pi^2 * h_bar^2).
    # The final value is the negative of this.
    lindhard_value = - (m_e * k_F / (pi**2 * h_bar**2))
    
    # Step 6: Print the components of the final equation and the result.
    # The equation is Π(0,0) = - (m_e * k_F) / (π^2 * ħ^2)
    print("\nThe Lindhard function at the specified limit is given by the equation:")
    print("Π(q=0, ω=0) = - (m_e * k_F) / (π² * ħ²)\n")
    print("Substituting the numerical values (in atomic units):")
    print(f"  m_e = {m_e:.1f}")
    print(f"  k_F = {k_F:.4f} (calculated from r_s)")
    print(f"  π   = {pi:.4f}")
    print(f"  ħ   = {h_bar:.1f}")
    
    print(f"\nΠ(q=0, ω=0) = - ({m_e:.1f} * {k_F:.4f}) / ({pi:.4f}² * {h_bar:.1f}²)")
    print(f"\nFinal numerical value (in atomic units): {lindhard_value:.4f}")

calculate_lindhard_function_zero_limit()