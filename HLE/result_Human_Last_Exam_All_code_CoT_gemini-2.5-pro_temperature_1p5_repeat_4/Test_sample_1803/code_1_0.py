def solve_energy_shift():
    """
    Calculates and prints the formula for the ground state zero-point energy
    shift for two interacting quantum harmonic oscillators.
    """
    # Symbolic representation of constants and variables
    e = 'e'  # electron charge
    hbar = 'hbar' # reduced Planck constant
    epsilon_0 = 'epsilon_0' # vacuum permittivity
    m = 'm' # mass
    omega_0 = 'omega_0' # angular frequency
    R = 'R' # distance

    # Coefficients and powers from the derivation
    num_coeff = -1
    e_pow = 4
    hbar_pow = 1
    den_coeff1 = 8
    four_pi_pow = 2
    m_pow = 2
    omega_0_pow = 3
    R_pow = 6
    
    # The number '4' is part of the 4*pi*epsilon_0 term
    four_in_paren = 4

    print("The leading term for the ground state energy shift (Delta_E) due to the Coulomb interaction is:")
    
    # Print the equation in a readable, single-line format
    equation = f"Delta_E = ({num_coeff} * {e}^{e_pow} * {hbar}) / ({den_coeff1} * ({four_in_paren}*pi*{epsilon_0})^{four_pi_pow} * {m}^{m_pow} * {omega_0}^{omega_0_pow} * {R}^{R_pow})"
    print("\n" + equation + "\n")

    print("Here is a breakdown of the numbers in the final equation:")
    print(f"  - The overall numerical coefficient in the numerator is: {num_coeff}")
    print(f"  - The power of the elementary charge 'e' is: {e_pow}")
    print(f"  - The power of the reduced Planck constant 'hbar' is: {hbar_pow}")
    print(f"  - The numerical coefficient in the denominator is: {den_coeff1}")
    print(f"  - The number within the 'permittivity' term is: {four_in_paren}")
    print(f"  - The power of the (4*pi*epsilon_0) term is: {four_pi_pow}")
    print(f"  - The power of the mass 'm' is: {m_pow}")
    print(f"  - The power of the angular frequency 'omega_0' is: {omega_0_pow}")
    print(f"  - The power of the distance 'R' is: {R_pow}")

if __name__ == '__main__':
    solve_energy_shift()