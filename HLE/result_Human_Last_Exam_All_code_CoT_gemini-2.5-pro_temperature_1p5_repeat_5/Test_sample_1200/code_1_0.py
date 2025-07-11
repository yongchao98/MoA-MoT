def analyze_alterations():
    """
    Models genetic alterations to demonstrate copy-neutral loss of heterozygosity.
    """
    # Start with a normal, heterozygous chromosome pair at a specific locus.
    # The cell has one paternal allele ('A') and one maternal allele ('B').
    # Copy number is 2. The cell is heterozygous.
    initial_chromosomes = ['Paternal_Allele_A', 'Maternal_Allele_B']

    print(f"Initial State: {initial_chromosomes}")
    print(f"Initial Copy Number: {len(initial_chromosomes)}")
    print("-" * 40)

    # --- Dictionary of Alterations and their simulated results ---
    alterations = {
        "A. Mitotic Recombination": ['Paternal_Allele_A', 'Paternal_Allele_A'],
        "B. A deletion of a chromosomal region": ['Paternal_Allele_A'],
        "C. Trisomy": ['Paternal_Allele_A', 'Maternal_Allele_B', 'Maternal_Allele_B'],
        "D. Uniparental Disomy": ['Maternal_Allele_B', 'Maternal_Allele_B'],
        "E. Duplication of a chromosomal region": ['Paternal_Allele_A', 'Maternal_Allele_B', 'Paternal_Allele_A']
    }

    # --- Analyze each alteration ---
    for name, resulting_alleles in alterations.items():
        print(f"Simulating: {name}")

        copy_number = len(resulting_alleles)
        # LOH occurs if only one type of allele remains (the set of unique alleles has a size of 1)
        is_loh = len(set(resulting_alleles)) == 1
        is_copy_neutral = copy_number == 2

        print(f"    Resulting Alleles: {resulting_alleles}")
        print(f"    Copy Number: {copy_number} ({'Neutral' if is_copy_neutral else 'Altered'})")
        print(f"    Loss of Heterozygosity (LOH): {is_loh}")
        
        # Check for the specific condition requested by the user
        if is_loh and is_copy_neutral:
            print("    CONCLUSION: This is Copy-Neutral Loss of Heterozygosity.")
        else:
            print("    CONCLUSION: This is NOT copy-neutral LOH.")
        print("-" * 40)

# Run the analysis
analyze_alterations()