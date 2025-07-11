def analyze_alteration(name, description, initial_locus, final_locus):
    """
    Analyzes a genetic locus to see if an alteration causes copy-neutral LOH.

    Args:
        name (str): The name of the alteration (e.g., "A. Mitotic Recombination").
        description (str): A brief explanation of the alteration's effect.
        initial_locus (list): The starting heterozygous locus.
        final_locus (list): The resulting locus after the alteration.
    """
    print(f"--- Analyzing: {name} ---")
    print(f"Description: {description}")
    print(f"Initial Locus (Heterozygous): {initial_locus}")
    print(f"Final Locus State: {final_locus}")

    # Condition 1: Is it copy-neutral?
    # Gene dosage is maintained if the number of alleles is the same as the start.
    is_copy_neutral = len(final_locus) == len(initial_locus)
    print(f"1. Is it Copy-Neutral? (Allele count: {len(initial_locus)} -> {len(final_locus)}): {is_copy_neutral}")

    # Condition 2: Is there Loss of Heterozygosity (LOH)?
    # LOH occurs if the locus is no longer heterozygous (i.e., all alleles are the same).
    is_loh = False
    if len(final_locus) > 0:
        # Check if all elements in the list are identical to the first one.
        is_loh = all(allele == final_locus[0] for allele in final_locus)
    print(f"2. Is there Loss of Heterozygosity? (Alleles were different, are now identical): {is_loh}")

    # Final Verdict
    if is_copy_neutral and is_loh:
        print("Verdict: This alteration perfectly matches the description of copy-neutral LOH.\n")
    else:
        print("Verdict: This alteration does NOT meet the criteria for copy-neutral LOH.\n")

def main():
    """
    Simulates various genetic alterations to identify the one causing copy-neutral LOH.
    """
    initial_locus = ['A', 'a']  # 'A' and 'a' are alleles on homologous chromosomes.

    print("Investigating which genetic alteration leads to copy-neutral loss of heterozygosity...\n")

    # A. Mitotic Recombination
    analyze_alteration(
        "A. Mitotic Recombination",
        "A process that can result in a daughter cell becoming homozygous for a locus while remaining diploid.",
        initial_locus,
        ['A', 'A'] # One possible outcome for a daughter cell.
    )

    # B. A deletion of a chromosomal region
    analyze_alteration(
        "B. A deletion of a chromosomal region",
        "One of the alleles is physically lost from the chromosome, reducing the gene copy number.",
        initial_locus,
        ['A'] # The chromosome segment with 'a' is deleted.
    )

    # C. Trisomy
    analyze_alteration(
        "C. Trisomy",
        "The cell gains an entire extra chromosome, increasing the gene copy number.",
        initial_locus,
        ['A', 'a', 'A'] # Gain of a chromosome with allele 'A'.
    )

    # D. Uniparental Disomy
    analyze_alteration(
        "D. Uniparental Disomy",
        "The cell has the normal number of chromosomes (two), but both are inherited from a single parent.",
        initial_locus,
        ['A', 'A'] # Both chromosomes now carry the 'A' allele from one parent.
    )

    # E. Duplication of a chromosomal region
    analyze_alteration(
        "E. Duplication of a chromosomal region",
        "A segment of a chromosome is duplicated, increasing the gene copy number.",
        initial_locus,
        ['A', 'a', 'A'] # The region with the 'A' allele is duplicated.
    )
    
    print("--- Final Conclusion ---")
    print("Both Mitotic Recombination (A) and Uniparental Disomy (D) result in a state that is copy-neutral with LOH.")
    print("However, Uniparental Disomy (UPD) is the term that specifically DESCRIBES the resulting state itself: having two copies of a chromosome derived from a single parent.")
    print("Mitotic recombination is a MECHANISM that can cause this. Therefore, UPD is the most direct and encompassing answer for the phenomenon described.")


if __name__ == "__main__":
    main()