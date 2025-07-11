def analyze_genetic_alterations():
    """
    Analyzes different genetic alterations to find the one causing
    copy-neutral loss of heterozygosity (CN-LOH).
    """

    options = {
        'A': 'Mitotic recombination',
        'B': 'A deletion of a chromosomal region',
        'C': 'Trisomy',
        'D': 'Uniparental disomy',
        'E': 'Duplication of a chromosomal region'
    }

    analysis = {
        'A': "Causes LOH. Can be copy-neutral. A crossover event during mitosis can result in daughter cells that are homozygous for a gene region while remaining diploid. This is a valid mechanism for CN-LOH.",
        'B': "Causes LOH but is NOT copy-neutral. Deleting a chromosomal region reduces the gene dosage from two copies to one.",
        'C': "Is NOT copy-neutral. Trisomy means having three copies of a chromosome instead of two, which increases the gene dosage.",
        'D': "Causes LOH and is copy-neutral. This occurs when a cell has two copies of a chromosome from one parent and none from the other. The total chromosome number is correct (diploid), but all heterozygosity on that chromosome is lost. This is a classic example of CN-LOH.",
        'E': "Is NOT copy-neutral. Duplicating a chromosomal region increases the gene dosage."
    }

    print("Analyzing the options for Copy-Neutral Loss of Heterozygosity (CN-LOH):")
    print("-" * 60)

    correct_answer = None
    for key, value in options.items():
        print(f"Option {key}: {value}")
        print(f"Analysis: {analysis[key]}\n")
        # Uniparental disomy is the textbook definition of whole-chromosome CN-LOH.
        # While mitotic recombination (A) can also cause this, uniparental disomy (D)
        # is the most direct and comprehensive answer describing this state.
        if key == 'D':
            correct_answer = key

    print("-" * 60)
    print(f"Conclusion: Both Mitotic Recombination (A) and Uniparental Disomy (D) can cause CN-LOH.")
    print("However, Uniparental Disomy is the term that best describes the state of having two chromosomes from one parent, which is the quintessential example of whole-chromosome copy-neutral LOH.")
    print(f"\nThe best answer is {correct_answer}: {options[correct_answer]}.")

analyze_genetic_alterations()