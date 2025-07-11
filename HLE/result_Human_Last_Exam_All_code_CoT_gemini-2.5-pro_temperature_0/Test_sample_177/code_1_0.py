def solve_immunology_question():
    """
    This function analyzes the options to determine which receptor would enable
    T cells to act as antigen-presenting cells (APCs) and prints the solution.
    """
    # Define the answer choices
    choices = {
        'A': 'CD86 - A co-stimulatory molecule.',
        'B': 'CD80 - A co-stimulatory molecule.',
        'C': 'MHC class I - Presents endogenous antigens, already on T cells.',
        'D': 'TIM-4 - A receptor for clearing apoptotic cells.',
        'E': 'MHC class II - The key receptor for presenting exogenous antigens to helper T cells.'
    }

    # Explanation
    explanation = """
    To function as an antigen-presenting cell (APC) in the way that dendritic cells or macrophages do, a cell must be able to present processed external antigens to helper T cells. This specific function is mediated by Major Histocompatibility Complex (MHC) class II molecules.

    - T cells naturally express MHC class I, but not MHC class II.
    - MHC class II is the essential receptor that allows professional APCs to activate helper T cells, which orchestrate the wider adaptive immune response.
    - By engineering a T cell to express MHC class II, it would gain the ability to present external antigens, thereby functioning as an APC.
    - Other options like CD80 and CD86 are co-stimulatory molecules that provide a second signal for T cell activation but do not present the antigen itself.
    """

    correct_answer_key = 'E'
    correct_answer_description = choices[correct_answer_key]

    print("--- Immunology Question Analysis ---")
    print(explanation)
    print("\n--- Conclusion ---")
    print(f"The correct receptor to enable T cells to act as antigen-presenting cells is:")
    print(f"Answer: {correct_answer_key}. {correct_answer_description}")

# Execute the function to get the answer
solve_immunology_question()