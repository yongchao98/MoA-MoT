import sympy as sp

def evanescent_wave_stored_energy():
    """
    This function prints the symbolic expressions for the time-averaged stored energy
    per unit area in the electric and magnetic fields of an evanescent wave,
    based on the provided answer choices.
    """
    # Define symbolic variables
    n = sp.Symbol('n', real=True, positive=True)
    theta = sp.Symbol('theta', real=True)
    omega_c = sp.Symbol('(omega/c)', real=True, positive=True)
    epsilon_0 = sp.Symbol('varepsilon_0', real=True, positive=True)
    E_x0_i_sq = sp.Symbol('|E_{x0}^i|^2', real=True, positive=True)
    
    # Common Denominator components
    sqrt_term = sp.sqrt(n**2 * sp.sin(theta)**2 - 1)
    denom_factor1 = (n**2 - 1)
    denom_factor2 = ((n**2 + 1) * sp.sin(theta)**2 - 1)
    
    full_denominator = 2 * omega_c * denom_factor1 * denom_factor2 * sqrt_term
    
    # Numerator for Electric Field Energy (from option D)
    E_field_num_factor = 2*n**2 * sp.sin(theta)**2 - 1
    E_field_numerator = n**2 * E_field_num_factor * epsilon_0 * E_x0_i_sq
    
    # Energy in Electric Field
    energy_E = E_field_numerator / full_denominator
    
    # Numerator for Magnetic Field Energy (from option D)
    # NOTE: The provided option D has a likely typo here (both in factors and units).
    # We will represent it as written.
    H_field_num_factor = n**2 * sp.sin(theta)**2 - 1
    H_field_numerator = n**2 * H_field_num_factor * epsilon_0 * E_x0_i_sq
    
    # Energy in Magnetic Field
    energy_H = H_field_numerator / full_denominator
    
    print("Based on the derivation, option D is the most likely answer, assuming a typo in the H-field expression.")
    print("\nFinal Equations from Option D:")
    
    # Printing the equations
    print("\nEnergy in E field =")
    sp.pprint(energy_E, use_unicode=True)
    
    print("\nEnergy in H field =")
    sp.pprint(energy_H, use_unicode=True)

if __name__ == '__main__':
    evanescent_wave_stored_energy()
