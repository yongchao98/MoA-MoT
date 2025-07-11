def solve_immunology_question():
    """
    This function analyzes the options to determine which receptor would enable
    a T cell to function as an antigen-presenting cell (APC).
    """
    options = {
        'A': 'CD86',
        'B': 'CD80',
        'C': 'MHC class I',
        'D': 'TIM-4',
        'E': 'MHC class II'
    }

    explanation = """
    T cells are not naturally professional antigen-presenting cells (APCs).
    A key function of APCs is to process external antigens and present them
    to helper T cells. This specific task is performed by MHC class II molecules.
    While T cells already have MHC class I (for presenting internal antigens) and
    can be induced to express co-stimulatory molecules like CD80/CD86, they
    lack MHC class II. Therefore, to engineer a T cell to act as an APC,
    it must be made to express MHC class II.
    """

    correct_answer_key = 'E'

    print("Analysis:")
    print(explanation)
    print("Final Answer:")
    print(f"The correct receptor is '{options[correct_answer_key]}'.")

solve_immunology_question()