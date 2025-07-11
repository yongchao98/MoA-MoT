def solve_immunology_question():
    """
    Analyzes the molecular requirements for a T cell to act as an
    Antigen-Presenting Cell (APC) and prints the reasoning.
    """

    question = "What type of receptor could be used to enable T cells to act as antigen-presenting cells?"

    options = {
        "A": "CD86",
        "B": "CD80",
        "C": "MHC class I",
        "D": "TIM-4",
        "E": "MHC class II"
    }

    correct_answer_key = "E"

    print("Analyzing the question: To make a T cell act like a professional Antigen-Presenting Cell (APC), we must give it the ability to present external antigens to helper T cells.")
    print("-" * 20)

    # Step-by-step reasoning
    print("Step 1: The key function of a professional APC is presenting processed external antigens on its surface.")
    print("Step 2: This presentation to helper T cells is specifically mediated by Major Histocompatibility Complex (MHC) class II molecules.")
    print("Step 3: T cells normally recognize antigens but do not present external ones, as they lack MHC class II.")
    print("Step 4: While CD80 and CD86 are important co-stimulatory molecules, and MHC class I presents internal antigens, neither provides the core function of presenting external antigens like an APC.")
    print("Step 5: Therefore, engineering a T cell to express MHC class II would grant it this crucial APC capability.")
    print("-" * 20)

    # Final answer
    final_answer_description = options[correct_answer_key]
    print(f"The correct receptor to engineer into T cells is: {final_answer_description}")
    print(f"This involves the MHC class II pathway, distinct from the MHC class I pathway.")

# Execute the function to provide the answer
solve_immunology_question()