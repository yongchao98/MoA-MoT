def solve_medical_question():
    """
    This function presents and solves a multiple-choice question
    regarding the mechanism of magnesium on blood pressure.
    """
    question = "By which mechanism can magnesium supplementation help lower blood pressure?"
    
    options = {
        'A': "Through direct vasodilation",
        'B': "By protecting elastic fibers from calcium deposition",
        'C': "By increasing white matter and gray matter in the brain",
        'D': "By stimulating an inflammatory response",
        'E': "It raises systemic blood calcium levels"
    }

    correct_answer = 'A'
    
    explanation = (
        "Magnesium helps lower blood pressure primarily through direct vasodilation. "
        "It functions as a natural calcium channel blocker. Calcium influx into vascular "
        "smooth muscle cells causes them to contract, leading to vasoconstriction and higher "
        "blood pressure. Magnesium competes with calcium and blocks its entry into these cells, "
        "which promotes relaxation of the blood vessels (vasodilation), thus lowering blood pressure."
    )

    print(f"Question: {question}\n")
    print("Options:")
    for key, value in options.items():
        print(f"{key}. {value}")
    
    print("\n--- Analysis ---")
    print(f"The correct answer is '{correct_answer}'.")
    print(f"\nExplanation: {explanation}")

# Execute the function to display the answer
solve_medical_question()