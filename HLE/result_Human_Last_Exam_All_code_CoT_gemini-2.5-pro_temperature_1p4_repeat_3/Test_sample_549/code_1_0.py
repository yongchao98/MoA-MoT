def evaluate_quantum_conductivity_correction():
    """
    This function explains the derivation of the quantum correction to conductivity
    (weak localization) in a 3D bulk semiconductor and prints the final formula.
    """
    print("### Evaluation of Quantum Correction to Conductivity in a Bulk Semiconductor ###")
    print("\n--- Plan ---")
    print("1. Identify the physical phenomenon: Weak Localization in a 3D system.")
    print("2. Describe the mechanism: Constructive interference of time-reversed electron paths enhances backscattering.")
    print("3. Outline the calculation: A momentum-space integral over the 'Cooperon' propagator with physical cutoffs.")
    print("4. Present the final formula for the conductivity correction, δσ.")

    print("\n--- Derivation Outline ---")
    print("The quantum correction to conductivity, δσ, arises from weak localization. It is calculated by considering the constructive quantum interference between an electron wave traversing a closed path and its time-reversed counterpart.")
    print("This enhances the probability of an electron returning to its origin, which slightly hinders transport and reduces overall conductivity.")
    print("\nThe calculation involves an integral in momentum (q) space. For a 3-dimensional system, the integral is performed over a spherical shell defined by two cutoffs:")
    print(" - The lower momentum cutoff is q_min = 1/L_φ, where L_φ is the phase coherence length.")
    print(" - The upper momentum cutoff is q_max = 1/l, where l is the electron mean free path.")

    print("\n--- Final Formula ---")
    print("Performing the integration with the standard prefactors from the Kubo formalism yields the weak localization correction to conductivity in 3D:")

    # The final equation is printed clearly below, showing all constants and variables.
    print("\nThe formula is:")
    print("  δσ = - (e² / (2 * π² * ħ)) * (1/l - 1/L_φ)")
    
    # Explaining the terms in the equation
    print("\nWhere:")
    print("  δσ  : The quantum correction to the conductivity.")
    print("  e   : The elementary charge.")
    print("  π   : The mathematical constant Pi (3.14159...).")
    print("  ħ   : The reduced Planck constant (h / 2π).")
    print("  l   : The electron elastic mean free path.")
    print("  L_φ : The electron phase coherence length.")
    print("\nNote: Since electrons undergo many elastic collisions before losing phase coherence, L_φ >> l. This makes the term (1/l - 1/L_φ) positive, and the overall correction δσ is negative, as expected.")


if __name__ == "__main__":
    evaluate_quantum_conductivity_correction()