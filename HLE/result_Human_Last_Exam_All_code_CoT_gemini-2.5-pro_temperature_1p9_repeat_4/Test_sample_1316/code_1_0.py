def explain_bridges_experiment():
    """
    Explains the genetic event leading to X0 males in Drosophila experiments.
    """

    # Step 1: Define the core observation.
    offspring_genotype = "X0"
    offspring_description = "An unexpected male with a single X chromosome and no Y chromosome."

    # Step 2: Explain the gametes required to produce this offspring.
    # Based on Bridges' work, the non-disjunction happens in the female.
    female_gamete = "nullo-X egg (contains no sex chromosome)"
    male_gamete = "normal X-bearing sperm"
    formation_process = f"An {offspring_genotype} zygote is formed when a {female_gamete} is fertilized by a {male_gamete}."

    # Step 3: Identify the cellular process responsible.
    cause_event = "Non-disjunction of X chromosomes during meiosis in the female."
    
    # Step 4: Detail why Meiosis I is the specific event.
    meiosis_I_explanation = """Non-disjunction during Meiosis I in the female causes the homologous pair of X chromosomes to fail to separate.
This single event produces two types of abnormal eggs:
1. Eggs with two X chromosomes (XX)
2. Eggs with no X chromosomes (nullo-X)"""

    full_context = """The power of Bridges' discovery was that this single event (Meiosis I non-disjunction) could explain not only the X0 males (from nullo-X eggs)
but also the exceptional XXY females he observed (from XX eggs fertilized by Y sperm).
This makes Meiosis I non-disjunction the most specific and complete explanation."""

    print("### Analysis of Bridges' Experiment ###")
    print(f"Observed Offspring Genotype: {offspring_genotype} ({offspring_description})")
    print("\nRequired Gametic Contribution:")
    print(formation_process)
    print("\nUnderlying Causal Event:")
    print(cause_event)
    print("\nDistinguishing the Meiotic Stage:")
    print(meiosis_I_explanation)
    print(full_context)
    print("\nConclusion: The evidence points to non-disjunction of the X chromosome in female meiosis I.")

explain_bridges_experiment()