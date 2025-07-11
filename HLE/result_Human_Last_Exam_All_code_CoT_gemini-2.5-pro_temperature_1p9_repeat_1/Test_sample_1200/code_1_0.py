import textwrap

def analyze_genetic_alterations():
    """
    Analyzes different genetic alterations to identify the one that causes
    copy-neutral loss of heterozygosity (LOH).
    """

    choices = {
        'A': {
            'term': 'Mitotic recombination',
            'explanation': 'This process can result in a daughter cell that is homozygous for a part of a chromosome while the overall chromosome number remains at two (diploid). This is a valid mechanism for copy-neutral LOH, but typically affects only a segment of the chromosome.',
            'is_correct': False
        },
        'B': {
            'term': 'A deletion of a chromosomal region',
            'explanation': 'This causes Loss of Heterozygosity (LOH) by removing one allele. However, it is NOT copy-neutral because the total number of gene copies is reduced from two to one. This is a copy number loss.',
            'is_correct': False
        },
        'C': {
            'term': 'Trisomy',
            'explanation': 'This is the gain of a whole chromosome, resulting in three copies instead of two. This is a copy number GAIN, not a copy-neutral event.',
            'is_correct': False
        },
        'D': {
            'term': 'Uniparental disomy',
            'explanation': 'This is the inheritance of two copies of a chromosome from a single parent. The cell remains diploid (two chromosomes total), so it is COPY-NEUTRAL. However, since both chromosomes originate from one parent, the cell becomes homozygous for all genes on that chromosome, which is a classic example of LOH. This mechanism perfectly fits the description of maintaining gene dosage while losing heterozygosity.',
            'is_correct': True
        },
        'E': {
            'term': 'Duplication of a chromosomal region',
            'explanation': 'This leads to an increase in the number of copies for a part of a chromosome. This is a copy number GAIN, not a copy-neutral event.',
            'is_correct': False
        }
    }

    print("Analysis of Genetic Alterations for Copy-Neutral Loss of Heterozygosity:\n")
    correct_answer_key = None
    for key, value in choices.items():
        print(f"Choice {key}: {value['term']}")
        # Wrap the explanation text for better readability
        wrapped_explanation = textwrap.fill(value['explanation'], width=80, initial_indent='  ', subsequent_indent='  ')
        print(wrapped_explanation)
        print("-" * 80)
        if value['is_correct']:
            correct_answer_key = key

    if correct_answer_key:
        print(f"\nConclusion: The best answer is {correct_answer_key}, {choices[correct_answer_key]['term']}.")
        print("It is the only option that describes an event that is, by definition, both copy-neutral and results in a loss of heterozygosity.")

analyze_genetic_alterations()
<<<D>>>