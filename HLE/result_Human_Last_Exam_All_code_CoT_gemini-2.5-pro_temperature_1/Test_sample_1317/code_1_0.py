def solve_recombination_frequency():
    """
    Determines the location of the highest recombination frequency in an
    interrupted mating experiment based on gene transfer order.
    """

    # 1. The gene order and direction of transfer are established from the problem.
    # The transfer starts at 'thr' and proceeds towards 'gal'.
    gene_transfer_order = ['thr', 'azi', 'gal']

    print("Step 1: Determine the order of gene transfer.")
    print(f"Based on the experiment, the order of transfer is: {gene_transfer_order[0]} -> {gene_transfer_order[1]} -> {gene_transfer_order[2]}")
    print("-" * 30)

    # 2. State the principle of interrupted mating experiments.
    print("Step 2: Apply the principle of interrupted mating.")
    print("In interrupted mating, genes are transferred linearly over time.")
    print("The frequency of recombination is highest for genes that are transferred earliest,")
    print("as they have the most time to be successfully transferred before the mating bridge breaks.")
    print("-" * 30)

    # 3. Identify the first gene transferred.
    first_gene = gene_transfer_order[0]
    print(f"Step 3: Identify the first gene transferred.")
    print(f"The first gene to enter the recipient cell is '{first_gene}'.")
    print("-" * 30)

    # 4. Conclude the location of the highest recombinant frequency.
    print("Step 4: Determine the location of highest recombinant frequency.")
    print(f"Therefore, the '{first_gene}' marker will appear in the highest number of recombinants.")
    print("This corresponds to the location on the chromosome immediately after the entry point of the thr+ marker.")
    print("-" * 30)
    
    print("\nFinal Conclusion: The highest frequency of recombinants should be expected 'Immediately after thr+'.")

solve_recombination_frequency()