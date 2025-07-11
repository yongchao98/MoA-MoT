import textwrap

def explain_cn_loh_question():
    """
    This function explains the reasoning behind the correct answer to the biology question.
    """
    print("Analysis of Genetic Alterations and Copy-Neutral Loss of Heterozygosity (CN-LOH)")
    print("=" * 70)

    # Define the core concept
    definition = (
        "The question asks to identify the genetic alteration that leads to "
        "copy-neutral loss of heterozygosity (CN-LOH). CN-LOH is an event where a cell "
        "that was heterozygous (e.g., has alleles 'A' and 'a') becomes homozygous "
        "(e.g., 'A' and 'A' or 'a' and 'a') for a gene or chromosome, without changing "
        "the total number of chromosomes. The gene dosage is maintained."
    )
    print(textwrap.fill(definition, width=70))
    print("\n--- Evaluating the Choices ---\n")

    # Analyze each choice
    choices = {
        'A': ('Mitotic recombination', 'This is a valid mechanism for CN-LOH. A crossover during mitosis can result in a daughter cell with two identical copies of a chromosome arm, which is a copy-neutral LOH event. It is a strong candidate.'),
        'B': ('A deletion of a chromosomal region', 'This causes LOH but is NOT copy-neutral. It reduces the copy number from two to one, which is a copy number loss.'),
        'C': ('Trisomy', 'This is the presence of three chromosomes instead of two. It is a copy number GAIN, not neutral.'),
        'D': ('Uniparental disomy (UPD)', 'This is when an individual has two homologous chromosomes from a single parent. By definition, this is a copy-neutral state (disomy). It results in LOH and is the classic example of CN-LOH.'),
        'E': ('Duplication of a chromosomal region', 'This results in three or more copies of genes in a region. It is a copy number GAIN, not neutral.')
    }

    for option, (title, explanation) in choices.items():
        print(f"{option}. {title}:")
        print(textwrap.fill(f"   - {explanation}", width=70, initial_indent='', subsequent_indent='     '))
        print()

    # Provide the conclusion
    conclusion = (
        "Conclusion: Both Mitotic Recombination (A) and Uniparental Disomy (D) are related to CN-LOH. "
        "However, Uniparental Disomy is the most encompassing term that describes the state of having two "
        "chromosomes from a single parental source, which perfectly fits the definition of copy-neutral loss of "
        "heterozygosity. Therefore, it is the best answer."
    )
    print("--- Conclusion ---\n")
    print(textwrap.fill(conclusion, width=70))


if __name__ == "__main__":
    explain_cn_loh_question()