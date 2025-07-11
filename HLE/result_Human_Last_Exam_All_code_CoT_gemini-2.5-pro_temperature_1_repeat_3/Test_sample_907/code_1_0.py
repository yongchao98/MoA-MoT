def display_absorption_equations():
    """
    This function prints the equations for the absorption cross-section
    for a chain of molecules under two different assumptions.
    """

    # --- Introduction and variable definitions ---
    print("This script provides the equations for the absorption cross-section σ(ωL).")
    print("The following variables are used:")
    print("  σ(ωL): Absorption cross-section as a function of laser frequency ωL")
    print("  C:     A constant containing fundamental physical constants")
    print("  N:     The number of molecules in the chain")
    print("  μ_eg:  The transition dipole moment of a single molecule")
    print("  ω_eg:  The electronic transition frequency of a single molecule")
    print("  ωL:    The carrier frequency of the laser pulse")
    print("  τ:     The duration of the Gaussian laser pulse")
    print("  E_eg:  The electronic transition energy of a single molecule (E_eg = ħω_eg)")
    print("  J:     The near-neighbor coupling energy")
    print("  ħ:     The reduced Planck's constant")
    print("-" * 50)

    # --- Case a) No interaction ---
    print("Case a) The interaction between molecules can be neglected.")
    print("The absorption spectrum is the sum of the spectra of N independent molecules.")
    # Equation using unicode for better readability
    equation_a = "σ_a(ωL) = C * N * |μ_eg|² * ωL * exp[-(ω_eg - ωL)² * τ² / 2]"
    print("\nThe equation for the absorption cross-section is:")
    print(equation_a)
    print("-" * 50)

    # --- Case b) Near-neighbor interaction ---
    print("Case b) The interaction between near-neighbors is considered (Frenkel Exciton Model).")
    print("The absorption is dominated by a single, collective 'bright' exciton state (k=0).")
    # Equation using unicode for better readability
    equation_b = "σ_b(ωL) = C * N * |μ_eg|² * ωL * exp[-((E_eg + 2*J)/ħ - ωL)² * τ² / 2]"
    print("\nThe equation for the absorption cross-section is:")
    print(equation_b)
    print("-" * 50)


# Execute the function to print the results.
display_absorption_equations()

# The final answer in the required format
final_answer = """
Case a) No interaction:
σ_a(ωL) = C * N * |μ_eg|² * ωL * exp[-(ω_eg - ωL)² * τ² / 2]

Case b) Near-neighbor interaction:
σ_b(ωL) = C * N * |μ_eg|² * ωL * exp[-((E_eg + 2*J)/ħ - ωL)² * τ² / 2]
"""
print(f"<<<{final_answer}>>>")