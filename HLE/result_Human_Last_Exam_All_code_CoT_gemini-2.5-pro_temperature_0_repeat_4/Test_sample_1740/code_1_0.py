def count_human_let7_family_members():
    """
    Calculates and displays the number of let-7 family members in humans,
    based on the count of distinct encoding genes.
    """
    # The let-7 family members are typically counted by their distinct precursor genes.
    # Some subfamilies (like let-7a and let-7f) are encoded by multiple genes.
    # miR-98 is also considered a member of the let-7 family.
    subfamily_gene_counts = {
        "let-7a": 3,
        "let-7b": 1,
        "let-7c": 1,
        "let-7d": 1,
        "let-7e": 1,
        "let-7f": 2,
        "let-7g": 1,
        "let-7i": 1,
        "miR-98": 1
    }

    print("The number of human let-7 family members is determined by counting the distinct genes that encode them.")
    print("The breakdown by sub-family is as follows:")
    for subfamily, count in subfamily_gene_counts.items():
        print(f"- Genes for {subfamily}: {count}")

    # Calculate the total and create the equation string
    total_members = sum(subfamily_gene_counts.values())
    equation_parts = [str(count) for count in subfamily_gene_counts.values()]
    equation_str = " + ".join(equation_parts)

    print("\nThe final calculation is:")
    print(f"{equation_str} = {total_members}")

    print(f"\nTo date, {total_members} members of the let-7 family, based on gene count, have been identified in humans.")

if __name__ == "__main__":
    count_human_let7_family_members()