def analyze_genetic_alterations():
    """
    This script analyzes different genetic alterations to determine which one causes
    copy-neutral loss of heterozygosity (CN-LOH) while maintaining gene dosage.
    """
    print("Objective: Identify the genetic alteration that leads to copy-neutral loss of heterozygosity (CN-LOH).\n")
    print("Let's break down the key terms:")
    print("1. Loss of Heterozygosity (LOH): A cell goes from having two different alleles (e.g., 'Aa') to having only one type ('AA' or 'aa').")
    print("2. Copy-Neutral: The total number of chromosomes remains two. Gene dosage is maintained.\n")

    print("--- Evaluating the Answer Choices ---\n")

    # Option A
    print("A. Mitotic recombination:")
    print("   - Analysis: This can cause LOH and is copy-neutral. However, it is a mechanism that can lead to Uniparental Disomy in a daughter cell. Let's see if a better option exists.\n")

    # Option B
    print("B. A deletion of a chromosomal region:")
    print("   - Analysis: This causes LOH but is NOT copy-neutral. The copy number is reduced from 2 to 1 for the deleted region. This is incorrect.\n")

    # Option C
    print("C. Trisomy:")
    print("   - Analysis: This is having three chromosomes, which is a copy number GAIN, not neutral. This is incorrect.\n")

    # Option D
    print("D. Uniparental disomy (UPD):")
    print("   - Analysis: This term means having two homologous chromosomes that both originated from the same parent.")
    print("   - Is it copy-neutral? Yes, there are two copies of the chromosome, so gene dosage is maintained.")
    print("   - Does it cause LOH? Yes, because both chromosomes are identical copies from one parent, any heterozygous loci become homozygous.")
    print("   - Conclusion: This perfectly matches the definition of copy-neutral loss of heterozygosity. This is the correct answer.\n")

    # Option E
    print("E. Duplication of a chromosomal region:")
    print("   - Analysis: This is a copy number GAIN for that specific region, not copy-neutral. This is incorrect.\n")

    print("--- Final Determination ---")
    print("The genetic alteration most likely to lead to copy-neutral loss of heterozygosity is Uniparental Disomy.")
    print("The correct choice is D.")

# Run the analysis
analyze_genetic_alterations()