def interrupted_mating_analysis():
    """
    This function explains the reasoning for determining the location of the
    highest recombination frequency in an E. coli interrupted mating experiment.
    """
    # The order of gene transfer indicates the position relative to the origin of transfer.
    gene_order = "thr -> azi -> gal"
    first_gene = "thr"
    
    print("--- Analysis of Bacterial Chromosome Mapping ---")
    print("\nPrinciple of Interrupted Mating:")
    print("In Hfr conjugation, genes are transferred linearly from a starting point (origin).")
    print("The mating is often interrupted, so genes closer to the origin are transferred more frequently.")
    print("Recombination frequency is highest for the genes that are transferred most frequently.")
    
    print("\nApplying the Principle to the Problem:")
    print(f"The observed order of gene transfer is: {gene_order}.")
    print(f"This means '{first_gene}' is the first gene to be transferred from the donor to the recipient.")
    print(f"Because '{first_gene}' is transferred first, it will be found in the largest number of recipient cells.")
    
    print("\nConclusion:")
    print(f"The highest frequency of recombination will occur at the locus of the first transferred gene, which is '{first_gene}'.")
    print("Therefore, the genetic location with the highest frequency of recombinants is at or immediately following the thr+ marker in the direction of transfer.")

interrupted_mating_analysis()