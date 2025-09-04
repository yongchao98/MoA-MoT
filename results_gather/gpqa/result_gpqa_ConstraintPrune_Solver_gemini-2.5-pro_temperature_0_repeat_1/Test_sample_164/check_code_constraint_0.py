def check_polymer_catalysis_answer():
    """
    Checks the correctness of the answer to the ethylene polymerization question.

    The function establishes a ground truth for each statement based on chemical facts
    and then evaluates the provided answer.
    """
    llm_answer = "C"
    question_context = "The core of the question is about the specific technology that enables the creation of *regular* branches from a single ethylene feedstock. This points to a highly selective oligomerization reaction running in tandem with polymerization."

    # Ground truth for each statement's validity and relevance
    statements_ground_truth = {
        'A': {
            'is_factually_correct': True,
            'relevance': "Partially relevant. Noble metals like Ni and Pd can oligomerize ethylene, but they are often less selective for a single oligomer compared to the best-in-class systems, making them less ideal for creating 'regular' branches. The statement is true but describes a suboptimal or different technology (e.g., SHOP process produces a distribution of olefins)."
        },
        'B': {
            'is_factually_correct': True,
            'relevance': "Relevant, but as a consequence, not a cause. This statement confirms the technology's success and implementation (e.g., by Chevron Phillips Chemical) but does not explain the underlying chemical principle that makes it possible."
        },
        'C': {
            'is_factually_correct': True,
            'relevance': "Most relevant. This statement points directly to the key chemical innovation. Chromium (Cr), a Group VIa metal, is the basis for catalysts that are exceptionally selective for the trimerization of ethylene to 1-hexene. This selectivity is crucial for producing the 'regular branches' mentioned in the question. This is the core enabling science."
        },
        'D': {
            'is_factually_correct': False,
            'relevance': "N/A. The statement is factually incorrect. Aluminum-based activators (co-catalysts) like alkylaluminums (e.g., TEA, TIBA) and methylaluminoxane (MAO) are critical for activating both the oligomerization (e.g., Cr-based) and polymerization catalysts in these systems."
        }
    }

    # 1. Check if the chosen answer is factually correct.
    if not statements_ground_truth[llm_answer]['is_factually_correct']:
        reason = f"Incorrect. The chosen answer '{llm_answer}' is factually false. "
        reason += statements_ground_truth[llm_answer]['relevance']
        return reason

    # 2. Among the factually correct statements (A, B, C), check if the chosen answer is the most relevant one.
    # The best answer in a scientific context is the one that explains the fundamental principle.
    is_most_relevant = statements_ground_truth[llm_answer]['relevance'].startswith("Most relevant")

    if is_most_relevant:
        return "Correct"
    else:
        # This handles the case where the LLM might have chosen A or B.
        reason = f"Incorrect. While statement '{llm_answer}' is factually true, it is not the best answer for the question asked. "
        reason += f"The chosen answer is only '{statements_ground_truth[llm_answer]['relevance'].lower()}' "
        reason += "Statement C is the most appropriate answer because it describes the core chemical technology that enables the selective process."
        return reason

# Execute the check
result = check_polymer_catalysis_answer()
print(result)