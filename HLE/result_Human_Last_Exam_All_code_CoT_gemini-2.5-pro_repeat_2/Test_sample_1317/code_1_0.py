def solve_mating_experiment():
    """
    Analyzes an interrupted mating experiment to find the location of highest recombinant frequency.
    """
    # Gene order on the chromosome
    gene_order = ["thr", "azi", "gal"]

    # The first gene to be transferred from the Hfr donor
    first_gene_transferred = "thr"

    print("Principle of Interrupted Mating:")
    print("In Hfr conjugation, genes are transferred linearly from an origin point.")
    print("The frequency of recombination for a gene is highest for the gene that is transferred first,")
    print("and the frequency decreases for genes that are transferred later.\n")

    print(f"Given Gene Order: {' - '.join(gene_order)}")
    print(f"First gene transferred is: {first_gene_transferred}\n")

    print("Reasoning:")
    print(f"1. Since '{first_gene_transferred}' is the first gene transferred, it is closest to the origin of transfer.")
    print(f"2. Therefore, '{first_gene_transferred}' will be successfully transferred to the most recipient cells before mating is interrupted.")
    print("3. The gene 'azi' will be transferred less frequently, and 'gal' will be transferred least frequently.")
    print("4. As recombination frequency is proportional to transfer frequency, the highest frequency of recombinants will be for the 'thr' gene.\n")

    print("Conclusion:")
    print("The location with the highest frequency of recombinants is the locus for the 'thr' gene.")
    print("This corresponds to choice A: Immediately after thr+.\n")

    final_answer = "A"
    print(f"<<<{final_answer}>>>")

solve_mating_experiment()