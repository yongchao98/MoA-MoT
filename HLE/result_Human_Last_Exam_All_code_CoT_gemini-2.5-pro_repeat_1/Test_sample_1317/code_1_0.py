def solve_recombination_frequency():
    """
    Explains the location of the highest recombination frequency in an
    interrupted mating experiment based on gene transfer order.
    """
    
    # Gene order based on time of entry
    gene_order = ["thr", "azi", "gal"]
    
    print("Principle of Interrupted Mating and Recombination:")
    print("In an Hfr x F- cross, genes are transferred from the donor to the recipient in a linear sequence.")
    print("The mating process is fragile and can be interrupted at any time.")
    print("Therefore, genes transferred earlier have a higher probability of being successfully integrated into the recipient's chromosome.")
    print("This results in a higher frequency of recombinants for earlier genes.\n")

    print("Experiment Details:")
    # The '->' represents the direction and order of transfer
    print(f"The observed gene transfer order is: {gene_order[0]} -> {gene_order[1]} -> {gene_order[2]}\n")

    # Determine the location with the highest frequency
    earliest_gene = gene_order[0]
    
    print("Conclusion:")
    print(f"The first gene to be transferred is '{earliest_gene}'.")
    print(f"Because '{earliest_gene}' is transferred most frequently, the region associated with this gene will show the highest frequency of recombinants.")
    print(f"This location is immediately after the '{earliest_gene}+' marker enters the recipient cell.")

solve_recombination_frequency()