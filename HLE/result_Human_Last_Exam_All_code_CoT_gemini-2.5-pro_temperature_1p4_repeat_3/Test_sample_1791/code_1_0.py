import sys

def solve_electron_phonon_interaction():
    """
    This function performs a symbolic derivation of the effective electron-electron 
    interaction mediated by phonons using the path integral method and prints the steps.
    """

    # --- Introduction and Plan ---
    print("### Plan to derive the effective electron-electron interaction ###")
    print("1. Start with the partition function in the path integral formalism, with the action separated into electron, phonon, and interaction parts.")
    print("2. Express the action in Matsubara frequencies.")
    print("3. The part of the action involving phonon fields is quadratic, which allows for exact integration (Gaussian integral).")
    print("4. 'Integrate out' the phonon fields by completing the square.")
    print("5. The result is a new term in the action that depends only on the electron density operators. This is the effective electron-electron interaction.")
    print("6. From this new term, we extract the effective interaction potential V_eff(q, ν_n).")
    print("7. Finally, we will print the resulting expression, as requested, for a single momentum mode 'q' without summing over it.")
    print("-" * 60)

    # --- Step-by-step Derivation ---
    print("\n### Step-by-step Derivation ###\n")

    # 1. The Action in Matsubara Space
    print("Step 1: The Action in Matsubara Frequency Space")
    print("The total action is S = S_el + S_ph + S_el-ph.")
    print("We focus on the phonon (S_ph) and electron-phonon interaction (S_el-ph) parts.")
    print("After Fourier transforming to Matsubara frequencies (τ -> ν_n), the action for the phonon fields (φ) is:")
    print("S_ph + S_el-ph = Σ_{q,j,n} [ φ_bar(q,j,n) * (-i*ν_n + w_q) * φ(q,j,n) + ...")
    print("                      ... + M(q,j) * ρ(-q,-n) * (φ(q,j,n) + φ_bar(-q,j,-n)) ]")
    print("where M(q,j) = g * i * q_j / sqrt(2 * m * w_q) is the electron-phonon coupling vertex.")
    print("ρ(q,n) is the electron density operator at momentum q and Matsubara frequency ν_n.")
    print("-" * 60)

    # 2. Rearranging the action
    print("\nStep 2: Rearranging the Action into a Standard Gaussian Integral Form")
    print("We rearrange the interaction term to group terms with φ and φ_bar.")
    print("The action can be written as: S_p = Σ_{q,j,n} [ φ_bar * A * φ + J_bar * φ + φ_bar * J ]")
    print("where:")
    print("  A = (-i*ν_n + w_q)")
    print("  J_bar(q,j,n) = M(-q,j) * ρ(q,n)")
    print("  J(q,j,n)     = M(q,j) * ρ(-q,-n)")
    print("-" * 60)

    # 3. Performing the Gaussian Integral
    print("\nStep 3: Integrating out the Phonon Fields")
    print("The integral ∫ d(φ_bar,φ) exp(-S_p) is a standard Gaussian integral.")
    print("The result of integrating over the phonon fields φ is exp(ΔS_eff), where ΔS_eff is the induced effective action for the electrons.")
    print("The formula for this integral gives: ΔS_eff = Σ_{q,j,n} J_bar * A⁻¹ * J")
    print("Substituting the terms from Step 2:")
    print("ΔS_eff = - Σ_{q,j,n} [ M(-q,j) * ρ(q,n) ] * [ 1 / (-i*ν_n + w_q) ] * [ M(q,j) * ρ(-q,-n) ]")
    print("Note: The minus sign comes from bringing the term to the exponent of the partition function e^(-S).")
    print("-" * 60)
    
    # 4. Finding the Effective Potential
    print("\nStep 4: Extracting the Effective Interaction Potential V_eff(q, ν_n)")
    print("The effective action has the form: ΔS_eff = (1/2) * Σ_{q,n} V_eff(q, ν_n) * ρ(q,n) * ρ(-q,-n)")
    print("By symmetrizing our result for ΔS_eff over (q, n) and (-q, -n), we find:")
    print("V_eff(q, ν_n) = - Σ_j M(q,j) * M(-q,j) * [ 1/(w_q - i*ν_n) + 1/(w_q + i*ν_n) ]")
    print("Simplifying the term in the square brackets gives: 2*w_q / (w_q² + ν_n²)")
    print("-" * 60)

    # 5. Final Result
    print("\nStep 5: Substituting the Coupling Constant to get the Final Result")
    print("Now we substitute the expressions for the coupling constants:")
    print("M(q,j) = g * i * q_j / sqrt(2*m*w_q)")
    print("M(-q,j) = g * i * (-q_j) / sqrt(2*m*w_q) = -M(q,j)")
    print("So, M(q,j) * M(-q,j) = -[M(q,j)]² = - (g² * (i*q_j)² / (2*m*w_q)) = g² * q_j² / (2*m*w_q)")
    
    print("\nSubstituting this into the expression for V_eff(q, ν_n):")
    print("V_eff(q, ν_n) = - Σ_j [ g² * q_j² / (2*m*w_q) ] * [ 2*w_q / (w_q² + ν_n²) ]")
    print("The terms '2*w_q' cancel out. The sum Σ_j q_j² is just the squared magnitude of the vector q, |q|².")
    print("-" * 60)
    
    # --- Final Expression ---
    print("\n### Final Result: The Effective Electron-Electron Interaction ###\n")
    print("The effective interaction potential between electrons, mediated by phonons, for a specific momentum transfer 'q' and energy transfer 'i*ν_n', is:")

    # Define symbolic variables for the final equation
    g = "g"
    q_vec_sq = "Σ_j(q_j)²" # which is also |q|²
    m = "m"
    w_q = "w_q"
    nu_n_sq = "ν_n²"
    rho_q = "ρ(q)"
    rho_minus_q = "ρ(-q)"

    # Construct and print the final equation for V_eff
    numerator = f"{g}² * {q_vec_sq}"
    denominator = f"{m} * ({w_q}² + {nu_n_sq})"
    print(f"               {numerator}")
    print(f"V_eff(q, ν_n) = - --------------------")
    print(f"               {denominator}")
    
    print("\nThis potential leads to an effective interaction term in the Hamiltonian for mode q:")
    
    # Construct and print the effective Hamiltonian term
    print(f"H_eff(q) = V_eff(q) * {rho_q} * {rho_minus_q}")
    
    print("\nNote: This interaction is frequency-dependent. At zero frequency (ν_n = 0), it is V_eff = -(g²|q|²)/(m*w_q²), which is always attractive.")
    print("This attraction is the mechanism behind conventional (BCS) superconductivity.")
    
# Execute the function
solve_electron_phonon_interaction()
# The final result is the formula derived and printed above.
# The structure of the answer is the python code itself.
# Let's construct the final output for the specified format.
# The core result is the expression for V_eff.
final_expression = "- (g**2 * Sum(q_j**2, j)) / (m * (w_q**2 + nu_n**2))"
# The question format seems to expect a single value or simple formula. 
# This derivation is the answer. I will put the main expression in the final tags.
final_answer = "-g^2 * |q|^2 / (m * (w_q^2 + v_n^2))"
sys.stdout = open(sys.stdout.fileno(), mode='w', encoding='utf-8', closefd=False) # allow unicode printing
# This print is for the final answer format, not for the user to see directly in the output block.
# print(f'<<<{final_answer}>>>')
