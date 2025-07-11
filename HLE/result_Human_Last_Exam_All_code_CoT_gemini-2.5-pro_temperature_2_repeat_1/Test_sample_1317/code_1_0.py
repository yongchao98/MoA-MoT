def explain_recombination_frequency():
    """
    Explains where the highest frequency of recombinants is expected
    in an E. coli interrupted mating experiment.
    """

    # The fundamental principle of interrupted mating experiments.
    principle = (
        "In E. coli conjugation, the chromosome is transferred linearly from a donor to a recipient.\n"
        "Genes transferred earlier enter a larger number of recipient cells than genes transferred later.\n"
        "The frequency of recombination for a gene is highest for the gene that is transferred most frequently."
    )

    # Information given in the problem.
    gene_order = "thr-azi-gal"
    first_marker = "thr+"

    print("--- Analysis of the Experiment ---")
    print(principle)

    print("\n--- Applying the Principle to the Problem ---")
    print(f"The established gene order is: {gene_order}")
    print(f"The first genetic marker observed to be transferred is: {first_marker}")

    # Step-by-step reasoning.
    print("\nReasoning:")
    print("1. Since 'thr' is the first gene to be transferred, it is closest to the origin of transfer.")
    print("2. Because it is transferred first, the 'thr' gene will be successfully received by the largest percentage of recipient cells.")
    print("3. As recombination frequency is proportional to transfer frequency, the highest frequency of recombinants will be for the 'thr' gene.")
    print("\nConclusion:")
    print("The highest frequency of recombinants should be expected at the genetic location of the first marker transferred, which is 'thr'.")
    print("This corresponds to the location immediately at and following the entry of the 'thr+' marker.")

explain_recombination_frequency()