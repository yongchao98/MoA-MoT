def check_alteration(name, alleles):
    """
    Analyzes a list of alleles to check for copy number and LOH.

    Args:
        name (str): The name of the genetic alteration.
        alleles (list): The list of alleles after the alteration.
    """
    copy_number = len(alleles)
    # LOH occurs if the number of unique alleles becomes 1.
    is_loh = len(set(alleles)) == 1
    is_copy_neutral = (copy_number == 2)

    print(f"--- Analysis for: {name} ---")
    print(f"Initial state: ['Allele_A', 'Allele_a']")
    print(f"Final state: {alleles}")
    print(f"Is Copy-Neutral (copy number = 2)? {is_copy_neutral}")
    print(f"Has Loss of Heterozygosity (LOH)? {is_loh}")
    if is_copy_neutral and is_loh:
        print("Result: This alteration causes Copy-Neutral LOH.")
    else:
        print("Result: This alteration does not meet the criteria for Copy-Neutral LOH.")
    print("-" * 30)

def main():
    """
    Simulates different genetic alterations to identify which one causes
    copy-neutral loss of heterozygosity (LOH).
    """
    # A. Mitotic Recombination can result in a cell that is homozygous.
    alleles_A = ['Allele_a', 'Allele_a']
    check_alteration("A. Mitotic Recombination", alleles_A)

    # B. Deletion of a chromosomal region removes an allele.
    alleles_B = ['Allele_A']
    check_alteration("B. Deletion", alleles_B)

    # C. Trisomy adds an extra chromosome.
    alleles_C = ['Allele_A', 'Allele_a', 'Allele_a']
    check_alteration("C. Trisomy", alleles_C)

    # D. Uniparental Disomy is the state of having two chromosomes from one
    #    parent. One mechanism is chromosome loss followed by duplication.
    #    e.g., lose 'Allele_a', then duplicate 'Allele_A'.
    alleles_D = ['Allele_A', 'Allele_A']
    check_alteration("D. Uniparental Disomy", alleles_D)

    # E. Duplication of a chromosomal region.
    alleles_E = ['Allele_A', 'Allele_a', 'Allele_a']
    check_alteration("E. Duplication", alleles_E)

    print("\nConclusion:")
    print("Both Mitotic Recombination (A) and Uniparental Disomy (D) result in Copy-Neutral LOH.")
    print("However, the question specifies 'maintaining the gene dosage despite allele deletion'.")
    print("This phrase best describes the mechanism of chromosome loss (an 'allele deletion' event) followed by duplication to restore gene dosage, which results in Uniparental Disomy.")
    print("Therefore, Uniparental Disomy is the most precise answer.")

if __name__ == "__main__":
    main()