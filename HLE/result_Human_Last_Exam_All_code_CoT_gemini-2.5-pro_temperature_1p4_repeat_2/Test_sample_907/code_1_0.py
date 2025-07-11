def solve_absorption_equation():
    """
    This function generates and prints the equations for the absorption
    cross-section for a chain of molecules under two different conditions.
    """

    # --- Introduction ---
    print("The absorption cross-section, sigma(omega), is calculated using first-order time-dependent perturbation theory.")
    print("The general form is a sum over all allowed transitions from an initial occupied state 'i' to a final unoccupied state 'f'.\n")
    print("-" * 50)

    # --- Symbolic Components of the Equations ---
    # These strings represent the physical quantities in the equations.
    sigma_omega = "sigma(omega)"
    proportionality_const = "C"
    frequency = "omega"
    num_molecules = "N"
    transition_dipole_if = "|mu_if|^2"
    lineshape_func = "G"
    hbar_omega = "hbar*omega"
    energy_diff_if = "(E_f - E_i)"
    coupling_const = "J_if"

    # --- Case a) Negligible Interaction ---
    print("\na) The interaction between molecules can be neglected.\n")
    print("In this case, the total absorption is the incoherent sum of the absorptions of N independent molecules.")
    print("The final equation is:")

    # Construct and print the equation for case (a)
    # The term Sum_if[...] represents the sum over all i->f transitions for a single molecule.
    # The total cross-section is N times this single-molecule contribution.
    equation_a = (
        f"{sigma_omega} = {proportionality_const} * {frequency} * {num_molecules} * "
        f"Sum_if[ {transition_dipole_if} * {lineshape_func}({hbar_omega} - {energy_diff_if}) ]"
    )
    print(f"\n{equation_a}\n")

    # --- Case b) Near-Neighbor Interaction ---
    print("-" * 50)
    print("\nb) The interaction between near-neighbors should be considered.\n")
    print("Interaction leads to the formation of delocalized excitons. An optical selection rule dictates that only the k=0 exciton is created.")
    print("This shifts the energy of each transition by an amount 2*J_if, where J_if is the coupling for that specific transition.")
    print("The final equation is:")

    # Construct and print the equation for case (b)
    # The structure is similar, but the energy term in the lineshape function is modified.
    energy_term_b = f"({energy_diff_if} + 2*{coupling_const})"
    equation_b = (
        f"{sigma_omega} = {proportionality_const} * {frequency} * {num_molecules} * "
        f"Sum_if[ {transition_dipole_if} * {lineshape_func}({hbar_omega} - {energy_term_b}) ]"
    )
    print(f"\n{equation_b}\n")

    # --- Explanation of Symbols ---
    print("-" * 50)
    print("\nExplanation of the terms:")
    print(f"  {sigma_omega:<15}: The absorption cross-section as a function of light frequency omega.")
    print(f"  {'C':<15}: A constant of proportionality containing fundamental constants (e.g., c, epsilon_0).")
    print(f"  {'N':<15}: The total number of molecules in the chain.")
    print(f"  {'Sum_if[...]':<15}: A sum over all initial occupied states 'i' and final unoccupied states 'f'.")
    print(f"  {transition_dipole_if:<15}: The squared transition dipole moment for a single molecule transition i -> f.")
    print(f"  {lineshape_func:<15}: A lineshape function (e.g., a Gaussian for a pulsed laser) describing the peak shape.")
    print(f"  {energy_diff_if:<15}: The transition energy for a single molecule.")
    print(f"  {coupling_const:<15}: The near-neighbor coupling energy for the specific i -> f transition (only in case b).")

solve_absorption_equation()