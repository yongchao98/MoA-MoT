def absorption_cross_section_equation():
    """
    This script explains and prints the equations for the absorption cross-section
    of a molecular chain under two different assumptions, based on first-order
    time-dependent perturbation theory for a Gaussian laser pulse.
    """
    
    print("=" * 80)
    print("Derivation of Absorption Cross-Section for a Molecular Chain")
    print("=" * 80)
    print("\nFrom first-order time-dependent perturbation theory, the absorption signal σ(ω)")
    print("for a system interacting with a Gaussian laser pulse is proportional to the material's")
    print("intrinsic absorption spectrum multiplied by the pulse's spectral shape.")
    
    print("\nThe general equation is:")
    print("σ(ω) ∝ ω * exp(-(ω - ω₀)²τ²/2) * Σ |<f|μ|G>|² * δ(ω - ω_f)")
    print("\nWhere:")
    print("  σ(ω)    : Absorption cross-section as a function of light frequency ω.")
    print("  ω       : Frequency of the absorbed light.")
    print("  ω₀      : Central frequency of the Gaussian laser pulse.")
    print("  τ       : Duration (standard deviation) of the Gaussian laser pulse.")
    print("  exp(...) : The Gaussian spectral shape of the laser pulse.")
    print("  |G>     : The ground state of the entire system.")
    print("  |f>     : A final (excited) state of the system.")
    print("  μ       : The total electric dipole moment operator.")
    print("  <f|μ|G> : The transition dipole moment between the ground and final state.")
    print("  ω_f     : The transition frequency, (E_f - E_G)/ħ.")
    print("  δ(...)  : The Dirac delta function, representing a sharp absorption line at ω = ω_f.")
    print("  Σ       : Sum over all possible final states |f>.")
    print("-" * 80)
    
    # --- Case a) No interaction between molecules ---
    print("\nCase a): The interaction between molecules can be neglected.\n")
    print("In this case, the molecules are independent. An excitation is localized on a single")
    print("molecule 'n'. There are N such possible excited states, |E_n>.")
    
    print("\n1. Energies:")
    print("All N excited states are degenerate. The transition frequency is the same for all,")
    print("equal to the single-molecule transition frequency, ω_mol.")
    print("  ω_f = ω_mol for all n=1..N")
    
    print("\n2. Transition Dipole Moments:")
    print("The total transition strength is the sum of the squared transition dipole moments (TDM)")
    print("to each of the N degenerate states. Since each state corresponds to exciting one")
    print("molecule, the TDM for each is just the TDM of a single molecule, μ_mol.")
    print("  Σ |<E_n|μ|G>|² = N * |μ_mol|²")
    
    print("\nFinal Equation for Case (a):")
    print("Combining these results, the sum collapses into a single term. The equation is:")
    print("\n  σ_a(ω) = C * N * |μ_mol|² * ω * exp(-(ω - ω₀)²τ²/2) * δ(ω - ω_mol)\n")
    print("Where:")
    print("  C        : Proportionality constant containing fundamental constants (ħ, c, ε₀).")
    print("  N        : The number of molecules in the chain.")
    print("  |μ_mol|² : The squared transition dipole moment of a single molecule.")
    print("  ω_mol    : The transition frequency of a single molecule, (E_elec_excited - E_elec_ground)/ħ.")
    print("-" * 80)
    
    # --- Case b) Near-neighbor interaction ---
    print("\nCase b): The interaction between near-neighbors should be considered.\n")
    print("Here, the electronic coupling J between adjacent molecules mixes the localized states |E_n>.")
    print("This creates new delocalized eigenstates called Frenkel excitons, |Ψ_k>.")
    
    print("\n1. Energies:")
    print("The degeneracy is lifted, and the N states form an 'exciton band'. For a ring model,")
    print("the energy of the exciton state |Ψ_k> is E_k = E_mol + 2J*cos(k), shifting the")
    print("transition frequencies.")
    print("  ω_k = ω_mol + 2J/ħ * cos(k)")
    
    print("\n2. Transition Dipole Moments (Selection Rule):")
    print("Due to symmetry, only the transition to the k=0 exciton state is allowed (it is 'bright').")
    print("For this state, all molecular dipoles add up coherently.")
    print("The transition dipole moment for this bright state is:")
    print("  |<Ψ_{k=0}|μ|G>|² = N * |μ_mol|²")
    print("For all other states (k≠0), the transition dipole moment is zero.")
    
    print("\nFinal Equation for Case (b):")
    print("The absorption spectrum is again a single peak, but its position is shifted by the coupling J.")
    print("\n  σ_b(ω) = C * N * |μ_mol|² * ω * exp(-(ω - ω₀)²τ²/2) * δ(ω - (ω_mol + 2J/ħ))\n")
    print("Where:")
    print("  J        : The near-neighbor coupling energy. If J < 0 (typical for H-aggregates),")
    print("             the absorption peak is blue-shifted. If J > 0 (J-aggregates), it is red-shifted.")
    print("  ħ        : The reduced Planck constant.")
    print("All other symbols are as defined in Case (a).")
    print("=" * 80)

# Execute the function to print the explanation
absorption_cross_section_equation()