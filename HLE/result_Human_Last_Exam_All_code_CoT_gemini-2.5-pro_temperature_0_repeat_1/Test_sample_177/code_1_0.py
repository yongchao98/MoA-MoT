def solve_immunology_question():
    """
    This function analyzes the options to determine which receptor would enable
    a T cell to act as an Antigen-Presenting Cell (APC).
    """
    options = {
        'A': 'CD86: A co-stimulatory molecule, provides the second signal for T cell activation but does not present the antigen.',
        'B': 'CD80: A co-stimulatory molecule, similar to CD86. It is crucial for activation but is not the antigen-presenting receptor.',
        'C': 'MHC class I: Presents endogenous (internal) antigens. T cells already express this. It is not used for presenting exogenous (external) antigens, which is the key role of a professional APC.',
        'D': 'TIM-4: A receptor for recognizing apoptotic cells, involved in phagocytosis, not antigen presentation to T cells.',
        'E': 'MHC class II: The key receptor used by professional APCs (like dendritic cells and macrophages) to present processed exogenous antigens to helper T cells. T cells do not normally express this.'
    }

    print("The goal is to engineer a T cell to act as an Antigen-Presenting Cell (APC).")
    print("The defining function of a professional APC is to present external antigens to helper T cells.")
    print("\nLet's analyze the options:\n")

    for key, description in options.items():
        print(f"Option {key}: {description}")

    correct_answer_key = 'E'
    print("\nConclusion:")
    print("To give a T cell the ability to present external antigens like a professional APC, it must be engineered to express the receptor responsible for this function.")
    print(f"Based on the analysis, this receptor is MHC class II. Therefore, the correct option is {correct_answer_key}.")

solve_immunology_question()