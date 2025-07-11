def demonstrate_multigene_family_compensation():
    """
    This script illustrates how multigene families, arising from gene duplication,
    can act as a compensatory mechanism for limited recombination.
    """
    # 1. Start with a single, essential gene in a hypothetical organism.
    # In a scenario with limited recombination, any deleterious mutation
    # to this gene could be permanently damaging to the lineage.
    genome = {'gene_101': 'vital_enzyme_production'}
    print("Initial State: A single essential gene.")
    print(f"Genome: {genome}\n")

    # 2. A gene duplication event occurs. This is a common type of mutation.
    # The result is a small "multigene family".
    genome['gene_101_copy_A'] = genome.pop('gene_101')
    genome['gene_101_copy_B'] = 'vital_enzyme_production' # The new copy
    print("Step 1: A gene duplication event occurs.")
    print(f"Genome now has a multigene family: {genome}\n")
    print("This provides redundancy. Both copies perform the same function.\n")

    # 3. One copy of the gene is now 'free' to accumulate mutations without
    # causing a loss of the original essential function, which is preserved by the other copy.
    # Let's say a series of mutations occurs in copy B.
    mutated_gene = 'gene_101_copy_B'
    original_function_gene = 'gene_101_copy_A'
    print(f"Step 2: Mutations accumulate in '{mutated_gene}'.")
    print(f"However, the essential function is safely preserved by '{original_function_gene}'.\n")

    # 4. Over evolutionary time, these mutations might result in a new, beneficial function.
    # This process is called neofunctionalization.
    genome[mutated_gene] = 'new_antibiotic_resistance'
    print("Step 3: The mutated gene evolves a new, beneficial function (neofunctionalization).")
    print(f"Final Genome: {genome}\n")

    print("Conclusion: The organism has not only retained its original vital function")
    print("but has also gained a new beneficial trait. This process of creating novel")
    print("genes via duplication and divergence compensates for the inability to create")
    print("new combinations of alleles through recombination.")

# Execute the demonstration
demonstrate_multigene_family_compensation()