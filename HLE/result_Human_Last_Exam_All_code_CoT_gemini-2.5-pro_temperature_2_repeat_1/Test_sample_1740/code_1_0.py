def count_human_let7_family_members():
    """
    Calculates and prints the total number of human let-7 family members,
    based on the number of distinct gene loci.

    The let-7 family includes several microRNAs, and some, like let-7a,
    are encoded by multiple distinct genes (paralogs). This function sums
    the gene counts for each family member to provide a total. miR-98 is
    included as it is widely considered a functional member of the let-7 family.
    """
    
    # A dictionary mapping each let-7 family member to the number of its
    # corresponding gene loci identified in the human genome.
    let7_gene_loci = {
        "let-7a": 3,
        "let-7b": 1,
        "let-7c": 1,
        "let-7d": 1,
        "let-7e": 1,
        "let-7f": 2,
        "let-7g": 1,
        "let-7i": 1,
        "miR-98": 1,
    }

    # Calculate the total number of members by summing the gene counts
    total_members = sum(let7_gene_loci.values())

    # Build the detailed equation string for clarity
    equation_parts = []
    for member, count in let7_gene_loci.items():
        equation_parts.append(f"{count} ({member})")
    
    equation_str = " + ".join(equation_parts)

    print(f"To date, the number of identified human let-7 family members (by gene loci) is calculated as follows:")
    # Print the full equation showing how the total is derived
    print(f"{equation_str} = {total_members}")


if __name__ == "__main__":
    count_human_let7_family_members()