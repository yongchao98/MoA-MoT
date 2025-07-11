def solve_tumorigenesis_question():
    """
    Analyzes genetic alterations to identify the one that causes
    copy-neutral loss of heterozygosity (CN-LOH).
    """
    # Step 1: Define the core concept of CN-LOH.
    # It has two main components:
    # 1. Loss of Heterozygosity (LOH): A cell that is heterozygous for a gene (e.g., has
    #    alleles 'A' and 'a') loses one allele and becomes homozygous (e.g., 'AA' or 'aa').
    # 2. Copy-Neutral: The total number of copies of the chromosome does not change.
    #    For a diploid cell, it remains at two copies.

    # Step 2: Define the answer choices and their characteristics.
    analysis = {
        'A': {
            'term': 'Mitotic recombination',
            'causes_LOH': True,
            'is_copy_neutral': True,
            'comment': 'This is a specific mechanism that can cause an exchange of genetic material, leading to homozygous daughter cells from a heterozygous parent without changing the chromosome count. This fits the definition of CN-LOH.'
        },
        'B': {
            'term': 'A deletion of a chromosomal region',
            'causes_LOH': True,
            'is_copy_neutral': False,
            'comment': 'This involves losing genetic material, which results in a copy number LOSS (from 2 copies to 1). It is therefore not copy-neutral.'
        },
        'C': {
            'term': 'Trisomy',
            'causes_LOH': False,
            'is_copy_neutral': False,
            'comment': 'This is the gain of a whole chromosome (from 2 to 3 copies). This is a copy number GAIN and is not copy-neutral.'
        },
        'D': {
            'term': 'Uniparental disomy',
            'causes_LOH': True,
            'is_copy_neutral': True,
            'comment': 'This describes the state where both homologous chromosomes come from a single parent. Acquired somatically (e.g., in a tumor), it means a cell loses one chromosome and duplicates the other, resulting in two identical copies. This is the very definition of CN-LOH.'
        },
        'E': {
            'term': 'Duplication of a chromosomal region',
            'causes_LOH': False,
            'is_copy_neutral': False,
            'comment': 'This involves repeating a segment of a chromosome, resulting in a copy number GAIN. It is not copy-neutral.'
        }
    }

    # Step 3 & 4: Evaluate choices and find the best answer.
    print("Evaluating Genetic Alterations for Copy-Neutral Loss of Heterozygosity (CN-LOH)\n")
    best_answer = None
    correct_choices = []
    
    for choice, data in analysis.items():
        if data['causes_LOH'] and data['is_copy_neutral']:
            correct_choices.append(choice)

    # In this case, both A and D are technically correct as they lead to CN-LOH.
    # However, Uniparental Disomy (D) is the broader and more definitive term.
    # Mitotic Recombination (A) is one specific *mechanism* that results in Uniparental Disomy.
    # Therefore, D is the most encompassing and accurate answer.
    best_answer = 'D'

    # Step 5: Print the detailed explanation.
    for choice, data in sorted(analysis.items()):
        print(f"--- Option {choice}: {data['term']} ---")
        print(f"Leads to Loss of Heterozygosity? {data['causes_LOH']}")
        print(f"Is Copy-Neutral? {data['is_copy_neutral']}")
        print(f"Conclusion: {data['comment']}\n")
    
    print("---------------------------------------------------------")
    print("Final Analysis:")
    print(f"The alterations that meet both criteria for CN-LOH are: {', '.join(correct_choices)} (Mitotic recombination and Uniparental disomy).")
    print("Uniparental disomy (D) is the most accurate and comprehensive answer because it is the direct description of the resulting state (two identical chromosomes), which can be caused by several mechanisms, including mitotic recombination (A).")
    print(f"\nThe most likely genetic alteration listed is D.")

solve_tumorigenesis_question()