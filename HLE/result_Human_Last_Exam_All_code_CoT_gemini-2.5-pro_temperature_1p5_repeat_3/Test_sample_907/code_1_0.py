def absorption_cross_section_equations():
    """
    This script prints the equations for the absorption cross-section
    of a molecular chain interacting with a Gaussian laser pulse,
    based on first-order time-dependent perturbation theory.
    """

    print("### Symbol Legend ###")
    print("---------------------------------------------------------------------")
    print("σ(ω)   : Absorption cross-section as a function of laser frequency ω")
    print("N      : Number of molecules in the chain")
    print("ω      : Central frequency of the Gaussian laser pulse (variable)")
    print("τ      : Duration of the Gaussian laser pulse")
    print("ħ      : Reduced Planck constant")
    print("c      : Speed of light")
    print("ε₀     : Vacuum permittivity")
    print("ε̂      : Polarization vector of the laser pulse")
    print("---------------------------------------------------------------------")
    print("For a single molecule:")
    print("ω_eg   : Transition frequency from ground state |g> to excited state |e>")
    print("μ_eg   : Transition dipole moment vector between |g> and |e>")
    print("---------------------------------------------------------------------")
    print("For the interacting chain:")
    print("J      : Near-neighbor coupling energy (a constant)")
    print("---------------------------------------------------------------------\n")


    # --- Case a) Non-interacting molecules ---
    print("### Case a) The interaction between molecules can be neglected ###\n")
    print("In this case, each molecule absorbs light independently. The total cross-section is simply N times the cross-section of a single molecule.")
    print("The absorption spectrum shows a single peak at the molecular resonance frequency ω_eg.\n")
    print("Final Equation:")
    print("σ_a(ω) = N * [ (2 * sqrt(π) * τ * ω_eg) / (ħ * c * ε₀) ] * |μ_eg ⋅ ε̂|² * exp[-(ω - ω_eg)² * τ²]\n")
    
    print("Breaking down the equation:")
    print("1. Overall Strength Factor: N")
    print("   - The absorption is proportional to the number of molecules.")
    print("2. Physical Constants and Pulse Parameters:")
    print("   - (2 * sqrt(π) * τ * ω_eg) / (ħ * c * ε₀)")
    print("3. Transition Strength (selection rule):")
    print("   - |μ_eg ⋅ ε̂|² (The squared projection of the transition dipole on the laser polarization)")
    print("4. Lineshape Function (from the Gaussian pulse):")
    print("   - exp[-(ω - ω_eg)² * τ²] (A Gaussian peak centered at ω_eg with width ~1/τ)")
    print("\n" + "="*70 + "\n")

    # --- Case b) Interacting molecules (Frenkel Exciton Model) ---
    print("### Case b) The interaction between near-neighbors should be considered ###\n")
    print("Interaction leads to collective, delocalized 'exciton' states. For a chain, selection rules dictate that only the k=0 exciton state is optically bright.")
    print("This state has its energy shifted by the coupling J, and its transition dipole is coherently enhanced.\n")

    # Define the new terms for case b
    k0_exciton_frequency = "ω_eg + 2*J/ħ"
    k0_exciton_dipole_sq = "N * |μ_eg ⋅ ε̂|²"
    
    print("Final Equation:")
    print(f"σ_b(ω) = [ (2 * sqrt(π) * τ * ({k0_exciton_frequency})) / (ħ * c * ε₀) ] * ({k0_exciton_dipole_sq}) * exp[-(ω - ({k0_exciton_frequency}))² * τ²]\n")

    print("Breaking down the equation (changes from Case a are highlighted):")
    print("1. Resonance Frequency (shifted):")
    print(f"   - The absorption peak is now at the k=0 exciton frequency: ω_k=0 = {k0_exciton_frequency}")
    print("   - This is a shift of 2*J/ħ from the single-molecule resonance. The sign of J determines if it's a blue-shift (J>0, H-aggregate) or red-shift (J<0, J-aggregate).")
    print("2. Transition Strength (coherently enhanced):")
    print(f"   - The squared transition dipole moment is now |μ_k=0 ⋅ ε̂|² = {k0_exciton_dipole_sq}")
    print("   - The intensity of the transition to the collective exciton state is N times that of a single molecule, a phenomenon known as superradiance.")
    print("3. Lineshape Function:")
    print(f"   - exp[-(ω - ({k0_exciton_frequency}))² * τ²] (A Gaussian peak centered at the new, shifted frequency)")
    print("\n" + "="*70)

if __name__ == '__main__':
    absorption_cross_section_equations()