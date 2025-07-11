def solve_history_question():
    """
    This function encapsulates the answer to the historical question
    based on the analysis provided.
    """
    
    # Components of the correct answer choice
    year_from_choice = 1223
    author_from_choice = "Rigord"
    final_answer_letter = 'C'

    # The prompt requests that numbers from any 'final equation' be printed.
    # As this is a historical question, we will interpret this by presenting
    # the numerical component of the chosen answer.
    print(f"The year in the selected answer is: {year_from_choice}")
    print(f"The author in the selected answer is: {author_from_choice}")

    # Present the final answer in the required format
    print(f"<<<{final_answer_letter}>>>")

solve_history_question()