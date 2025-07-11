def check_conditions(alleles, scenario_name):
    """Checks a list of alleles for LOH and copy-neutral status."""
    print(f"--- Evaluating: {scenario_name} ---")
    
    initial_alleles = ['WT', 'mutant']
    initial_copy_number = len(initial_alleles)
    
    # Condition 1: Loss of Heterozygosity (LOH)
    # This means the 'WT' allele is gone, and only one type of allele remains.
    is_loh = 'WT' not in alleles and len(set(alleles)) == 1
    
    # Condition 2: Copy-Neutral
    # The final number of alleles is the same as the initial number.
    is_copy_neutral = len(alleles) == initial_copy_number
    
    print(f"Initial state: {initial_alleles}, Copy Number: {initial_copy_number}")
    print(f"Final state:   {alleles}, Copy Number: {len(alleles)}")
    print(f"Is there Loss of Heterozygosity (LOH)? -> {is_loh}")
    print(f"Is the event Copy-Neutral? -> {is_copy_neutral}")
    
    if is_loh and is_copy_neutral:
        print("Result: This event leads to Copy-Neutral Loss of Heterozygosity.\n")
    else:
        print("Result: This event does NOT lead to Copy-Neutral Loss of Heterozygosity.\n")

def main():
    # We start with a heterozygous cell for a tumor suppressor gene.
    # ['WT', 'mutant'] -> WT = Wild Type, mutant = non-functional allele
    
    # A. Mitotic Recombination
    # A crossover event can lead to a daughter cell that inherits two 'mutant' chromatids.
    alleles_A = ['mutant', 'mutant']
    check_conditions(alleles_A, "A. Mitotic Recombination")

    # B. A deletion of a chromosomal region
    # The region containing the 'WT' allele is deleted.
    alleles_B = ['mutant']
    check_conditions(alleles_B, "B. A deletion of a chromosomal region")

    # C. Trisomy
    # The cell gains an extra chromosome, in this case one carrying the 'mutant' allele.
    alleles_C = ['WT', 'mutant', 'mutant']
    check_conditions(alleles_C, "C. Trisomy")

    # D. Uniparental Disomy
    # The chromosome with the 'WT' allele is lost, and the 'mutant' chromosome is duplicated.
    # The cell ends up with two copies of the 'mutant' chromosome.
    alleles_D = ['mutant', 'mutant']
    check_conditions(alleles_D, "D. Uniparental Disomy")
    
    # E. Duplication of a chromosomal region
    # The region containing the 'mutant' allele is duplicated.
    alleles_E = ['WT', 'mutant', 'mutant']
    check_conditions(alleles_E, "E. Duplication of a chromosomal region")

main()