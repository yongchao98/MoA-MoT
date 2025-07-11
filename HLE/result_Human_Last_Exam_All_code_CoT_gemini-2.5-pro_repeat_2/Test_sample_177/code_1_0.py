def solve_immunology_question():
    """
    This function analyzes the biological roles of different cellular receptors
    to determine which one would enable a T cell to act as an Antigen-Presenting Cell (APC).
    """

    question = "What type of receptor could be used to enable T cells to act as antigen-presenting cells?"

    options = {
        'A': 'CD86',
        'B': 'CD80',
        'C': 'MHC class I',
        'D': 'TIM-4',
        'E': 'MHC class II'
    }

    # Explanation of why MHC class II is the correct answer
    explanation = [
        "1. Professional Antigen-Presenting Cells (APCs) present external antigens to helper T cells.",
        "2. The specific molecule used for this presentation is the MHC class II receptor.",
        "3. T cells do not naturally express MHC class II and thus cannot act as professional APCs.",
        "4. Engineering a T cell to express MHC class II would grant it this specific capability.",
        "5. Other options are incorrect: CD80/CD86 are for co-stimulation, MHC class I is for internal antigens, and TIM-4 is for clearing apoptotic cells."
    ]

    correct_answer_key = 'E'
    correct_answer_value = options[correct_answer_key]

    print(f"Question: {question}\n")
    print("Analysis:")
    for step in explanation:
        print(f"- {step}")

    print("\n---")
    print(f"Final Answer: The correct choice is {correct_answer_key}, which is {correct_answer_value}.")
    print("---")


solve_immunology_question()