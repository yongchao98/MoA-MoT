def solve_immunology_question():
    """
    Analyzes the options to determine which receptor would enable T cells
    to act as Antigen-Presenting Cells (APCs).
    """
    choices = {
        'A': 'CD86',
        'B': 'CD80',
        'C': 'MHC class I',
        'D': 'TIM-4',
        'E': 'MHC class II'
    }

    correct_answer_key = 'E'
    correct_answer_value = choices[correct_answer_key]

    explanation = (
        "To enable a T cell to function as an Antigen-Presenting Cell (APC), it must be able to present "
        "exogenous (external) antigens to other T cells, specifically CD4+ helper T cells. This is the "
        "defining role of professional APCs.\n\n"
        "The molecular machinery for this function is the Major Histocompatibility Complex (MHC) class II receptor. "
        "While T cells naturally express MHC class I to present internal antigens, they do not express MHC class II.\n\n"
        "By engineering a T cell to express MHC class II, it would gain the ability to process and present "
        "external antigens, effectively allowing it to act as an APC.\n\n"
        f"Therefore, the correct choice is E: {correct_answer_value}."
    )

    print(explanation)

solve_immunology_question()