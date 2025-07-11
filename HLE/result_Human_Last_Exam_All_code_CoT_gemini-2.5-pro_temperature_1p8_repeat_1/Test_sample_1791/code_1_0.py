import sympy as sp

def derive_effective_interaction():
    """
    This function presents a step-by-step symbolic derivation of the effective 
    electron-electron interaction mediated by phonons. The derivation is based on 
    integrating out the phonon degrees of freedom from the total Hamiltonian 
    using a path integral formalism.
    """
    
    print("Derivation of the Effective Electron-Electron Interaction\n")

    # --- 1. Define symbolic variables ---
    print("Step 1: Define the symbolic variables.")
    # Physical constants and variables
    w_q = sp.Symbol("ω_q", real=True, positive=True) # Phonon frequency
    q_j = sp.Symbol("q_j", real=True)               # Phonon wave vector component
    g = sp.Symbol("g", real=True)                   # Electron-phonon coupling constant
    m = sp.Symbol("m", real=True, positive=True)   # Mass (e.g., electron mass)
    
    # Matsubara frequency for bosons
    w_n = sp.Symbol("ω_n", real=True)              
    
    # Operators in the Hamiltonian
    rho_q = sp.Symbol("ρ_q")                        # Electron density operator for wavevector q
    rho_mq_n = sp.Symbol("ρ_{-q}(-ω_n)")             # Fourier component of rho_{-q}
    rho_q_n = sp.Symbol("ρ_q(ω_n)")                  # Fourier component of rho_q
    a_q, a_mq_dag = sp.symbols("a_q, a_{-q}^†")      # Phonon operators

    print(f"  Phonon frequency: {w_q}")
    print(f"  Matsubara frequency: {w_n}")
    print(f"  Electron-phonon coupling constant: {g}")
    print(f"  Phonon momentum component: {q_j}")
    print(f"  Mass: {m}")
    print("-" * 50)

    # --- 2. Define the Interaction Hamiltonian ---
    print("Step 2: State the Electron-Phonon Interaction Hamiltonian.")
    # Define the coupling constant M_q from the problem statement
    M_q = g * sp.I * q_j / sp.sqrt(2 * m * w_q)
    print("The electron-phonon interaction Hamiltonian couples electron density ρ_q to the phonon field.")
    print("For a single mode q (as requested), the relevant part of the interaction is:")
    print(f"  H_el-ph(q) = M_q * {rho_q} * ({a_q} + {a_mq_dag})")
    print("  where the coupling strength M_q is given as:")
    display_M_q = sp.Eq(sp.Symbol("M_q"), M_q)
    print(f"  {display_M_q}")
    print("-" * 50)

    # --- 3. Outline the Path Integral procedure ---
    print("Step 3: Integrate out phonon fields.")
    print("To find the effective electron-electron interaction, we integrate out the phonon")
    print("degrees of freedom in the path integral formulation of the partition function.")
    print("This procedure is a Gaussian integration over the phonon fields, which are coupled")
    print("to the electron density acting as a source term.")
    print("The result is a new term in the effective action for electrons:")
    print(f"  S_eff = - (1/2) * Σ_{{q,ω_n}} |M_q|² * D(q, iω_n) * {rho_q_n} * {rho_mq_n}")
    print("where D(q, iω_n) is the propagator for the phonon field operator (a_q + a_{-q}^†).")
    print("-" * 50)

    # --- 4. State the Phonon Propagator ---
    print("Step 4: Define the Phonon Propagator D(q, iω_n).")
    # The propagator for the field phi_q = a_q + a_{-q}^dagger is given by:
    D_qw = 2 * w_q / (w_n**2 + w_q**2)
    print("The propagator for the phonon field (a_q + a_{-q}^†) in Matsubara frequency space is:")
    display_D_qw = sp.Eq(sp.Symbol("D(q, iω_n)"), D_qw)
    print(f"  {display_D_qw}")
    print("-" * 50)

    # --- 5. Calculate the effective interaction potential ---
    print("Step 5: Compute the effective interaction potential V_eff(q, iω_n).")
    print("The effective potential is the kernel of the interaction term: V_eff = -|M_q|² * D(q, iω_n).")

    # Calculate |M_q|^2
    M_q_sq_abs = (g**2 * (sp.I * q_j) * (-sp.I * q_j)) / (sp.sqrt(2 * m * w_q)**2)
    M_q_sq = sp.simplify(M_q_sq_abs)

    print("\nFirst, we calculate the magnitude squared of the coupling constant, |M_q|²:")
    display_M_q_sq = sp.Eq(sp.Abs(sp.Symbol("M_q"))**2, M_q_sq)
    print(f"  {display_M_q_sq}")
    
    # Now, calculate V_eff
    V_eff_qw = -M_q_sq * D_qw
    V_eff_qw = sp.simplify(V_eff_qw)
    
    print("\nNow, we multiply -|M_q|² by the propagator D(q, iω_n):")
    print(f"  V_eff = (-|M_q|²) * D(q, iω_n)")
    display_V_eff = sp.Eq(sp.Symbol("V_eff(q, iω_n)"), V_eff_qw)
    print(f"  Result: {display_V_eff}")
    print("-" * 50)
    
    # --- 6. Final Result: The Effective Interaction Term ---
    print("Step 6: State the final expression for the effective interaction.")
    print("The resulting phonon-mediated interaction between electrons, for a given mode q and")
    print("Matsubara frequency ω_n, is described by the potential V_eff(q, iω_n).")
    
    print("\n========================= FINAL EQUATION =========================")
    print("The effective electron-electron interaction potential V_eff is:")
    # Using pretty print for a more readable console output
    print(sp.pretty(display_V_eff, use_unicode=True))

    print("\nTo explicitly show the components of the final equation:")
    # Numerator parts
    num_part1 = -g**2
    num_part2 = q_j**2
    # Denominator parts
    den_part1 = m
    den_part2_w_n = w_n**2
    den_part2_w_q = w_q**2
    
    print(f"  V_eff = ( {num_part1} * {num_part2} ) / ( {den_part1} * ({den_part2_w_n} + {den_part2_w_q}) )")
    print("================================================================")


if __name__ == '__main__':
    # Execute the function to print the derivation
    derive_effective_interaction()