def analyze_genetic_alterations():
    """
    Analyzes different genetic alterations to determine which leads to
    copy-neutral loss of heterozygosity (LOH).
    """
    print("Analyzing the options for copy-neutral loss of heterozygosity (LOH):")

    print("\nKey Concepts:")
    print(" - Loss of Heterozygosity (LOH): A cell goes from being heterozygous (e.g., alleles A/a) to homozygous (e.g., A/A or a/a).")
    print(" - Copy-Neutral: The total number of chromosomes remains the same (e.g., stays diploid).")

    options = {
        'A': "Mitotic recombination: Can create a homozygous daughter cell from a heterozygous parent cell while maintaining a diploid state. This is a valid mechanism for copy-neutral LOH.",
        'B': "A deletion of a chromosomal region: Leads to LOH but is NOT copy-neutral. It results in a loss of genetic material and reduces the gene copy number from 2 to 1.",
        'C': "Trisomy: Involves gaining an extra chromosome. This is a copy number GAIN, not neutral.",
        'D': "Uniparental disomy (UPD): A cell has two copies of a chromosome from a single parent. This is by definition copy-neutral. When it results from duplicating one chromosome and losing the other (isodisomy), it causes LOH for the entire chromosome. This is a classic example of copy-neutral LOH.",
        'E': "Duplication of a chromosomal region: Involves gaining extra genetic material. This is a copy number GAIN, not neutral."
    }

    print("\nEvaluating each choice:")
    for choice, explanation in options.items():
        print(f"  {choice}. {explanation}")

    print("\nConclusion:")
    print("Both Mitotic Recombination (A) and Uniparental Disomy (D) can cause copy-neutral LOH.")
    print("However, Uniparental Disomy is the term that best describes the resulting state of having two chromosomes from one parent, which perfectly fits the definition of being copy-neutral while causing a loss of heterozygosity. It is the most precise answer.")

analyze_genetic_alterations()
print("\n<<<D>>>")