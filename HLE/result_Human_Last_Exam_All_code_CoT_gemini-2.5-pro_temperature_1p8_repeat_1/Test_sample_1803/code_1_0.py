import sympy

def solve_energy_shift():
    """
    Calculates the ground state energy shift for two interacting
    quantum harmonic oscillators using second-order perturbation theory.
    """
    # Define symbolic variables for the physical quantities.
    # We use 'e2' to represent e^2 as given in the potential form e^2/(4*pi*r).
    m, omega0, R, hbar = sympy.symbols('m omega_0 R hbar', positive=True, real=True)
    e2 = sympy.Symbol('e^2', positive=True, real=True)
    x1, x2 = sympy.symbols('x_1 x_2', real=True)

    # --- Step 1: Derive the leading term of the interaction potential ---
    # The full potential is V = e^2 * [1/R - 1/(R+x2) - 1/(R-x1) + 1/(R+x2-x1)].
    # Expanding this for R >> x1, x2, the leading non-zero term is the
    # dipole-dipole interaction potential.
    # (Detailed expansion is V ~ (e^2/R) * (-2*x1*x2/R^2) )
    V_interaction = -2 * e2 * x1 * x2 / R**3
    
    print("Step 1: Determine the Interaction Potential (Perturbation)")
    print("----------------------------------------------------------")
    print("By Taylor expanding the Coulomb potential for large distance R, the leading")
    print("term of the interaction potential V is found to be:")
    print(f"V = {sympy.pretty(V_interaction)}")
    print("\n")

    # --- Step 2: First-Order Perturbation ---
    print("Step 2: First-Order Energy Correction")
    print("--------------------------------------")
    print("The first-order energy correction is the expectation value of V in the ground state |00>.")
    print("ΔE^(1) = <00| V |00> = <0|x1|0> * <0|x2|0> * (coefficient)")
    print("Since the expectation value of position x for any QHO eigenstate is zero, <n|x|n> = 0.")
    print("Therefore, the first-order energy correction is ΔE^(1) = 0.")
    print("\n")

    # --- Step 3: Second-Order Perturbation ---
    print("Step 3: Second-Order Energy Correction")
    print("---------------------------------------")
    print("We calculate the second-order correction using the formula:")
    print("ΔE^(2) = Σ_{n,m≠0} |<nm|V|00>|^2 / (E_0 - E_nm)")
    print("\nThe matrix element <nm|V|00> is non-zero only for n=1 and m=1.")
    
    # The matrix element <1|x|0> for a QHO
    x_10 = sympy.sqrt(hbar / (2 * m * omega0))
    print(f"\nThe required matrix element <1|x|0> for a single oscillator is: {sympy.pretty(x_10)}")

    # The full matrix element <11|V|00>
    V_coeff = -2 * e2 / R**3
    V_11_00 = V_coeff * x_10 * x_10
    print(f"\nThe matrix element <11|V|00> is: {sympy.pretty(V_11_00.simplify())}")
    
    # The square of the matrix element
    V_11_00_sq = V_11_00**2
    print(f"\nThe square of this matrix element is: {sympy.pretty(V_11_00_sq.simplify())}")

    # The energy denominator
    # E_ground = E_00 = hbar*omega0
    # E_excited = E_11 = 3*hbar*omega0
    energy_denom = (hbar * omega0) - 3 * (hbar * omega0)
    print(f"\nThe energy denominator E_0 - E_11 is: {sympy.pretty(energy_denom)}")
    print("\n")
    
    # --- Step 4: Final Calculation ---
    print("Step 4: Final Energy Shift Calculation")
    print("--------------------------------------")
    
    # Calculate the final energy shift
    delta_E_2 = V_11_00_sq / energy_denom
    
    # Extract coefficients for clear output
    v_11_00_num_coeff = V_coeff.args[0] * V_coeff.args[1] # -2
    v_11_00_sq_num_coeff = (v_11_00_num_coeff * x_10.args[1].args[0] / x_10.args[1].args[1])**2 # (-2 * 1/2)**2 = (-1)**2 = 1
    
    e_den_num_coeff = energy_denom.args[0] # -2
    
    print("The final energy shift is given by the formula:")
    print("ΔE_0 ≈ ΔE^(2) = |<11|V|00>|^2 / (E_0 - E_11)")
    print("\nSubstituting the numerical coefficients:")
    print(f"Coefficient of |<11|V|00>|^2 term = ({int(v_11_00_num_coeff)} * (1/2))^2 = {int(v_11_00_sq_num_coeff)}")
    print(f"Coefficient of (E_0 - E_11) term = {int(e_den_num_coeff)}")
    
    final_coeff_num = v_11_00_sq_num_coeff
    final_coeff_den = e_den_num_coeff
    
    print(f"\nResulting numerical factor in the final equation = {int(final_coeff_num)} / {int(final_coeff_den)} = {final_coeff_num/final_coeff_den}")
    
    print("\nThe final expression for the ground state energy shift is:")
    final_expression = sympy.simplify(delta_E_2)
    # Use pretty print for better formatting of the final answer
    sympy.pprint(final_expression, use_unicode=True)
    
    # To conform to the final answer format
    return final_expression


if __name__ == '__main__':
    final_answer = solve_energy_shift()
    # The final answer is symbolic, so we present it as a string
    answer_string = str(final_answer)
    # The output format is <<<answer content>>>. 
    # For a complex symbolic answer, presenting the full string is the most accurate representation.
    # print(f"\n<<<{answer_string}>>>") 
    # The final print above gives a symbolic expression, which seems more useful than just a number.
    # However, let me try to adhere to the final format. Since the question is about calculation,
    # the process and the final formula is the answer. So I will put the final formula inside the desired format.
    final_pretty_string = sympy.sstr(final_answer, order=None)
    print(f"\n<<<{final_pretty_string}>>>")
