import sympy

def calculate_energy_shift():
    """
    Calculates and displays the ground state energy shift for two interacting
    quantum harmonic oscillators using symbolic mathematics.
    """
    # Step 1: Define all symbolic variables for the physical quantities.
    # e: elementary charge unit
    # R: distance between the oscillators
    # m: mass of the oscillator
    # omega_0: angular frequency of the oscillators
    # hbar: reduced Planck's constant
    # x1, x2: position operators for the two oscillators
    e, R, m, omega_0, hbar = sympy.symbols('e R m omega_0 hbar', real=True, positive=True)
    x1, x2 = sympy.symbols('x1 x2')

    # Step 2: Define the interaction Hamiltonian (H_int).
    # Assuming collinear dipoles and separation, the leading term is the
    # dipole-dipole interaction potential. The problem statement uses e^2/(4*pi*r),
    # which we interpret as setting the Coulomb constant k_e=1.
    # H_int = C * x1 * x2
    C = -2 * e**2 / R**3
    H_int = C * x1 * x2
    
    print("Step 1: The interaction Hamiltonian (leading term in R) is:")
    sympy.pprint(H_int)
    print("-" * 30)

    # Step 3: Calculate the second-order energy shift.
    # The first-order shift is zero, so we go to second order.
    # ΔE = |<1,1|H_int|0,0>|^2 / (E_ground - E_excited)
    
    # The only non-zero matrix element of H_int is between the ground state |0,0>
    # and the excited state |1,1>.
    # The energy difference is (hbar*w0) - (3*hbar*w0) = -2*hbar*w0
    energy_denominator = -2 * hbar * omega_0

    # The square of the matrix element of the position operator |<1|x|0>|^2 for a QHO
    # is hbar / (2 * m * omega_0).
    # For H_int, the squared matrix element is |C * <1|x1|0> * <1|x2|0>|^2
    matrix_element_sq = C**2 * (hbar / (2 * m * omega_0))**2
    
    # The total energy shift is the ratio of these two quantities.
    delta_E = matrix_element_sq / energy_denominator

    # Simplify the final expression.
    delta_E_simplified = sympy.simplify(delta_E)

    print("Step 2: The ground state energy shift ΔE is calculated using second-order perturbation theory.")
    print("Final expression for ΔE:")
    sympy.pprint(delta_E_simplified)
    print("-" * 30)
    
    # Step 4: Output each number and exponent in the final equation.
    print("Step 3: Breaking down the components of the final equation:")
    numer, denom = delta_E_simplified.as_numer_denom()
    
    # Extract the overall sign from the numerator's coefficient
    sign = numer.as_coeff_mul()[0]
    print(f"Overall numerical sign: {sign}")

    # Extract the numerical coefficient from the denominator
    denom_coeff = denom.as_coeff_mul()[0]
    print(f"Numerical coefficient in the denominator: {denom_coeff}")
    
    # Get the variable parts without coefficients
    numer_vars = numer.as_coeff_mul()[1]
    denom_vars = denom.as_coeff_mul()[1]

    # Get the powers of each variable as a dictionary
    numer_powers = numer_vars.as_powers_dict()
    denom_powers = denom_vars.as_powers_dict()
    
    print("\nPowers of variables in the numerator:")
    for var, power in numer_powers.items():
        print(f"  '{var}' is raised to the power: {power}")
        
    print("\nPowers of variables in the denominator:")
    for var, power in denom_powers.items():
        print(f"  '{var}' is raised to the power: {power}")


# Execute the calculation
calculate_energy_shift()