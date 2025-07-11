def display_absorption_equations():
    """
    This script prints the equations for the absorption cross-section
    for a chain of molecules under two different conditions.
    """

    print("This script provides the equations for the absorption cross-section (σ) for a molecular chain interacting with an ultrashort Gaussian laser pulse, based on first-order time-dependent perturbation theory.")
    print("-" * 80)

    # --- Case (a): Negligible interaction ---
    print("a) The interaction between molecules can be neglected.\n")
    print("In this case, the total absorption is the sum of independent absorptions from each molecule.")
    print("The absorption cross-section σ as a function of the laser's angular frequency ω is:\n")
    
    # The number '2' is explicitly part of the equation's denominator.
    equation_a = "σ(ω) = C_a * (E_f - E_i) * |d_fi|^2 * exp[ -((E_f - E_i) - ħω)^2 * τ^2 / (2 * ħ^2) ]"
    print(equation_a)
    
    print("\nWhere:")
    print("  - C_a: A proportionality constant.")
    print("  - E_i: Energy of the initial molecular state (ground state).")
    print("  - E_f: Energy of the final molecular state (excited state).")
    print("  - d_fi: The transition dipole moment between states |i> and |f>.")
    print("  - ħ: The reduced Planck's constant.")
    print("  - ω: The angular frequency of the laser.")
    print("  - τ: The duration of the Gaussian laser pulse.")
    print("  - The number 2 appears in the denominator of the exponent, characteristic of a Gaussian function.")
    
    print("\n" + "=" * 80 + "\n")

    # --- Case (b): Near-neighbor interaction ---
    print("b) The interaction between near-neighbors should be considered.\n")
    print("In this case, near-neighbor coupling (J) leads to the formation of delocalized Frenkel excitons.")
    print("Due to selection rules in a 1D chain, only the k=0 exciton state is optically active.")
    print("The absorption cross-section σ is given by:\n")
    
    # The number '2' appears as a factor of J and in the denominator.
    equation_b = "σ(ω) = C_b * N * ((E_f - E_i) + 2*J) * |d_fi|^2 * exp[ -(((E_f - E_i) + 2*J) - ħω)^2 * τ^2 / (2 * ħ^2) ]"
    print(equation_b)

    print("\nWhere:")
    print("  - C_b: A proportionality constant.")
    print("  - N: The number of molecules in the chain.")
    print("  - J: The near-neighbor coupling energy.")
    print("  - (E_f - E_i), d_fi, ħ, ω, τ are the same as defined in case (a).")
    print("  - The energy of the transition is shifted by 2*J.")
    print("  - The absorption strength is enhanced by a factor of N (coherent effect).")
    print("  - The number 2 appears as a factor multiplying the coupling energy J and in the denominator of the exponent.")

if __name__ == '__main__':
    display_absorption_equations()