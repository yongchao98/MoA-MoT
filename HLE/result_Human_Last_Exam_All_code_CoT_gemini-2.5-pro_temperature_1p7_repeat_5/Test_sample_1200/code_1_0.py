def analyze_loh_mechanisms():
    """
    Analyzes different genetic alterations to determine which causes
    copy-neutral loss of heterozygosity (LOH).
    """
    initial_genotype = "['A', 'b']"
    initial_copy_number = 2

    print("Analyzing which mechanism leads to copy-neutral Loss of Heterozygosity (LOH).\n")
    print(f"Initial State: A heterozygous locus with genotype {initial_genotype} and copy number {initial_copy_number}.\n")
    print("--------------------------------------------------")

    options = {
        'A': 'Mitotic recombination',
        'B': 'A deletion of a chromosomal region',
        'C': 'Trisomy',
        'D': 'Uniparental disomy',
        'E': 'Duplication of a chromosomal region'
    }

    # Analysis for each option
    # A. Mitotic Recombination
    print(f"Choice A: {options['A']}")
    print("   - Result: Can result in a daughter cell that is homozygous for the locus, e.g., ['A', 'A'].")
    print(f"   - Copy Number Change: {initial_copy_number} -> 2")
    print("   - Is it LOH? Yes.")
    print("   - Is it Copy-Neutral? Yes. This is a possible mechanism.\n")

    # B. Deletion
    print(f"Choice B: {options['B']}")
    print("   - Result: One allele is lost, resulting in ['A'] or ['b'].")
    print(f"   - Copy Number Change: {initial_copy_number} -> 1")
    print("   - Is it LOH? Yes (hemizygous).")
    print("   - Is it Copy-Neutral? No. Copy number is reduced.\n")

    # C. Trisomy
    print(f"Choice C: {options['C']}")
    print("   - Result: An extra chromosome is gained, e.g., ['A', 'A', 'b'].")
    print(f"   - Copy Number Change: {initial_copy_number} -> 3")
    print("   - Is it LOH? No, heterozygosity can be maintained.")
    print("   - Is it Copy-Neutral? No. Copy number is increased.\n")

    # D. Uniparental Disomy (UPD)
    print(f"Choice D: {options['D']}")
    print("   - Result: One parental chromosome is lost and the remaining one is duplicated. Genotype becomes ['A', 'A'].")
    print(f"   - Copy Number Change: {initial_copy_number} -> 2")
    print("   - Is it LOH? Yes.")
    print("   - Is it Copy-Neutral? Yes. This is a classic mechanism for copy-neutral LOH.\n")

    # E. Duplication
    print(f"Choice E: {options['E']}")
    print("   - Result: A region is duplicated, leading to a state like ['A', 'A', 'b'].")
    print(f"   - Copy Number Change: {initial_copy_number} -> 3")
    print("   - Is it LOH? No, heterozygosity is maintained.")
    print("   - Is it Copy-Neutral? No. Copy number is increased.\n")

    print("--------------------------------------------------")
    print("Conclusion:")
    print("Both Mitotic Recombination (A) and Uniparental Disomy (D) result in a copy number of 2 and can cause LOH.")
    print("However, Uniparental Disomy is the most direct and accurate term describing the specific event where an individual has two homologous chromosomes from one parent, perfectly matching the definition of copy-neutral LOH.")

analyze_loh_mechanisms()
<<<D>>>