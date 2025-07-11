def display_absorption_equations():
    """
    This script prints the equations for the absorption cross-section
    of a molecular chain interacting with a Gaussian laser pulse,
    derived from first-order time-dependent perturbation theory.
    """
    
    print("="*65)
    print("General Symbols Used in the Equations")
    print("="*65)
    print("  σ(ħω): Absorption cross-section as a function of photon energy ħω.")
    print("  C_const: Proportionality constant containing fundamental physical values.")
    print("  ω: Angular frequency of the laser pulse.")
    print("  τ: Duration (half-width at 1/e) of the Gaussian laser pulse.")
    print("  ħ: Reduced Planck constant.")
    print("  d: Electric dipole moment operator.")
    print("  ε: Polarization vector of the laser's electric field.")
    print("  exp[...]: The exponential function, representing the Gaussian lineshape.")
    print("\n")


    # --- Case (a): No interaction between molecules ---
    print("="*65)
    print("Case a): No Interaction Between Molecules")
    print("="*65)
    print("The system is treated as a collection of independent molecules. The electronic")
    print("states are localized on each molecule.")
    print("\n--- Equation for Absorption Cross-Section (Case a) ---")
    
    # Equation components for case (a)
    sum_term_a = "Σ_{i,f}"
    dipole_term_a = "|<φ_f| d·ε |φ_i>|²"
    energy_term_a = "(E_f - E_i - ħω)"
    gaussian_term_a = f"exp[-{energy_term_a}² * τ² / (2ħ²)]"
    
    print(f"\n  σ(ħω)  =  C_const * ω * {sum_term_a} {dipole_term_a} * {gaussian_term_a}\n")

    print("--- Explanation of Terms (Case a) ---")
    print(f"1. {sum_term_a}: Sum over all occupied initial molecular orbitals |φ_i> (e.g., HOMO)\n"
          "   and all unoccupied final molecular orbitals |φ_f> (e.g., LUMO).")
    print(f"2. {dipole_term_a}: The squared transition dipole moment. This term contains the\n"
          "   selection rules; it is non-zero only for allowed transitions.")
    print(f"3. {gaussian_term_a}: The Gaussian lineshape function from the pulse.")
    print(f"   - {energy_term_a}: The energy difference between the final state (E_f),\n"
          "     the initial state (E_i), and the incident photon (ħω).")
    print("   - The absorption is maximized when ħω ≈ E_f - E_i. The width of the\n"
          "     absorption peak is determined by τ (shorter pulse = wider peak).")
    print("\n")


    # --- Case (b): Nearest-neighbor interaction ---
    print("="*65)
    print("Case b): Nearest-Neighbor Interaction Considered")
    print("="*65)
    print("Interaction between molecules delocalizes the electronic states into energy")
    print("bands, described by a tight-binding model.")
    print("\n--- Equation for Absorption Cross-Section (Case b) ---")

    # Equation components for case (b)
    sum_term_b = "Σ_k"
    dipole_term_b = "|d_cv·ε|²"
    delta_e_k = "ΔE(k)"
    energy_term_b = f"({delta_e_k} - ħω)"
    gaussian_term_b = f"exp[-{energy_term_b}² * τ² / (2ħ²)]"

    print(f"\n  σ(ħω)  =  C_const * ω * {dipole_term_b} * {sum_term_b} {gaussian_term_b}\n")
    
    print("--- Explanation of Terms (Case b) ---")
    print(f"1. {dipole_term_b}: The single-molecule transition dipole moment. Due to the\n"
          "   k-selection rule (k_final = k_initial), this becomes a constant factor.")
    print(f"2. {sum_term_b}: Sum over all allowed wavevectors 'k' in the first Brillouin zone.")
    print(f"3. {gaussian_term_b}: The Gaussian lineshape for a k-dependent transition.")
    print(f"   - {delta_e_k}: The transition energy now depends on the wavevector 'k'.")
    print(f"   - The formula for this energy is:")
    print(f"       {delta_e_k} = ΔE_mol - 2*J_eff*cos(k*a)")
    print(f"         - ΔE_mol: Energy gap of an isolated molecule (e.g., E_LUMO - E_HOMO).")
    print(f"         - J_eff: Effective interaction energy between adjacent molecules.")
    print(f"         - a: Spacing between molecules in the chain.")
    print("   - The absorption is no longer a single peak but a band, reflecting the\n"
          "     density of states for transitions across the energy gap ΔE(k).")
    print("\n")

if __name__ == '__main__':
    display_absorption_equations()