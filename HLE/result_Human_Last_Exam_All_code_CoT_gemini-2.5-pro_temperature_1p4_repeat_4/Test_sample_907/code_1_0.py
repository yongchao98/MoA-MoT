def print_absorption_cross_section_equations():
    """
    This function prints the formatted equations for the absorption cross-section
    for a molecular chain under two different assumptions:
    a) No interaction between molecules.
    b) Nearest-neighbor interactions are included (Frenkel exciton model).
    """

    # --- Common Symbols ---
    omega_0 = "ω₀"  # Laser central frequency
    tau = "τ"       # Pulse duration
    hbar = "ħ"      # Reduced Planck's constant
    c = "c"         # Speed of light
    epsilon_0 = "ε₀" # Vacuum permittivity
    pi = "π"        # Pi
    N = "N"         # Number of molecules
    mu_eg = "μ_eg"    # Transition dipole moment of a single molecule
    omega_eg = "ω_eg" # Molecule transition frequency

    # --- Prefactor C is common to both equations ---
    # C = (2 * ω₀ * τ * sqrt(π)) / (ħ * c * ε₀)
    prefactor_C = f"(2 * {omega_0} * {tau} * sqrt({pi})) / ({hbar} * {c} * {epsilon_0})"

    # --- Case a) Non-interacting molecules ---
    print("---")
    print("Case a) Equation for absorption cross-section (non-interacting molecules):")
    print("---\n")

    # Equation: σ_a(ω₀) = C * N * |μ_eg|² * exp(-(ω_eg - ω₀)² * τ²)
    term_a = f"{N} * |{mu_eg}|²"
    exp_a = f"exp(-({omega_eg} - {omega_0})² * {tau}²)"
    equation_a = f"σ_a({omega_0}) = {prefactor_C} * {term_a} * {exp_a}"
    
    print(equation_a)
    print("\nThis equation describes a single absorption peak with a Gaussian shape.")
    print("\n" + "="*80 + "\n")


    # --- Case b) Near-neighbor interaction ---
    print("---")
    print("Case b) Equation for absorption cross-section (with near-neighbor interaction):")
    print("---\n")

    # Symbols specific to case b
    J = "J" # Near-neighbor coupling energy
    j = "j" # Exciton quantum number
    n = "n" # Molecule index in the chain

    # Equation: σ_b(ω₀) = C * Σ_{j=1..N} |d_{G,j}|² * exp(-(ω_{G,j} - ω₀)²τ²)
    sum_term_b = f"Σ_{{{j}=1..{N}}}"
    dipole_term_b = f"|d_{{G,{j}}}|²"
    exp_b = f"exp(-({omega_eg} + (2{J}/{hbar})*cos({j}{pi}/({N}+1)) - {omega_0})² * {tau}²)"
    
    # Simplified main equation line for clarity
    main_equation_b = f"σ_b({omega_0}) = {prefactor_C} * {sum_term_b} [ {dipole_term_b} * exp(-(ω_{{G,{j}}} - {omega_0})² * {tau}²) ]"

    print(main_equation_b)

    print("\nwhere the exciton transition frequency ω_{G,j} and transition dipole moment d_{G,j} are:")
    
    # Definition of ω_{G,j}
    omega_Gj_def = f"ω_{{G,{j}}} = {omega_eg} + (2{J}/{hbar}) * cos({j}{pi} / ({N}+1))"
    print(f"\n1. Exciton Transition Frequency:\n   {omega_Gj_def}")

    # Definition of |d_{G,j}|²
    dipole_Gj_def = f"|d_{{G,{j}}}|² = |{mu_eg}|² * (2/({N}+1)) * [ Σ_{{{n}=1..{N}}} sin({j}{pi}{n} / ({N}+1)) ]²"
    print(f"\n2. Squared Transition Dipole Moment (for excitons):\n   {dipole_Gj_def}")

    print("\nThis equation describes a spectrum of multiple absorption peaks, corresponding to different exciton states.")
    print("\n" + "="*80 + "\n")

    # --- Legend for all symbols ---
    print("Legend of Symbols:")
    print(f"  σ({omega_0})  : Absorption cross-section as a function of laser frequency")
    print(f"  {N}        : Number of molecules in the chain")
    print(f"  {J}        : Near-neighbor coupling energy (interaction strength)")
    print(f"  {mu_eg}    : Transition dipole moment of a single, isolated molecule")
    print(f"  {omega_eg}   : Electronic transition frequency of a single molecule")
    print(f"  {omega_0}    : Central frequency of the Gaussian laser pulse")
    print(f"  {tau}      : Temporal duration (width) of the Gaussian laser pulse")
    print(f"  {hbar}       : Reduced Planck's constant")
    print(f"  {c}        : Speed of light in vacuum")
    print(f"  {epsilon_0}  : Permittivity of free space")
    print(f"  {pi}       : The mathematical constant pi")

if __name__ == "__main__":
    print_absorption_cross_section_equations()
