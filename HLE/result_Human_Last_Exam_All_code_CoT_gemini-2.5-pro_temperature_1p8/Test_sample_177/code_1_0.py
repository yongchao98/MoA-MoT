def find_best_receptor_for_apc_function():
    """
    Analyzes potential receptors to engineer a T cell into an Antigen-Presenting Cell (APC).

    The model assigns a score to each receptor based on its primary function in
    the context of presenting exogenous antigens, which is the key role of an APC.
    """
    # Options provided in the problem
    options = {
        'A': 'CD86',
        'B': 'CD80',
        'C': 'MHC class I',
        'D': 'TIM-4',
        'E': 'MHC class II'
    }

    # Assigning relevance scores based on biological function for APC-like activity
    # 10: Primary molecule for presenting exogenous antigens
    # 7: Essential co-stimulatory molecule, but not for presentation itself
    # 4: Presents endogenous antigens, not the target function
    # 2: Involved in uptake, a prerequisite, but not presentation
    relevance_scores = {
        'CD86': 7,
        'CD80': 7,
        'MHC class I': 4,
        'TIM-4': 2,
        'MHC class II': 10
    }

    print("Evaluating receptors based on their function for enabling T cells to act as APCs.")
    print("The 'equation' involves finding the maximum relevance score among the options.")
    print("Here are the numbers (scores) used in our evaluation:")

    # Print each score, fulfilling the requirement to output each number
    for option_letter, receptor_name in options.items():
        score = relevance_scores[receptor_name]
        print(f"  - Option {option_letter} ({receptor_name}): Score = {score}")

    # Find the receptor with the highest score
    best_receptor = max(relevance_scores, key=relevance_scores.get)
    max_score = relevance_scores[best_receptor]

    # Find the corresponding letter for the best receptor
    final_answer_letter = [letter for letter, name in options.items() if name == best_receptor][0]

    print(f"\nThe equation to find the solution is: max({list(relevance_scores.values())}) = {max_score}")
    print(f"The highest score of {max_score} belongs to '{best_receptor}'.")
    print(f"\nConclusion: To enable a T cell to present exogenous antigens like a professional APC, it should be engineered to express {best_receptor}.")
    print(f"This corresponds to answer choice {final_answer_letter}.")

# Run the analysis
find_best_receptor_for_apc_function()