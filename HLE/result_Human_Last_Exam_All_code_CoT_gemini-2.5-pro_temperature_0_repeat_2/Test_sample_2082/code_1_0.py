def solve_medical_question():
    """
    This function analyzes the provided multiple-choice question about magnesium's effect on blood pressure
    and prints the correct answer with a detailed explanation.
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
        "Magnesium helps lower blood pressure primarily by acting as a natural calcium channel blocker. "
        "It competes with calcium for entry into vascular smooth muscle cells. The influx of calcium "
        "causes these muscle cells to contract, leading to vasoconstriction and higher blood pressure. "
        "By blocking calcium's entry, magnesium promotes the relaxation of these cells, resulting in "
        "vasodilation (the widening of blood vessels). This vasodilation reduces peripheral resistance, "
        "which in turn lowers blood pressure."
    )
    
    print(f"Question: {question}\n")
    print("Answer Choices:")
    for key, value in options.items():
        print(f"{key}. {value}")
    
    print("\n----------------------------------")
    print(f"Selected Answer: {correct_answer}")
    print(f"Mechanism: {options[correct_answer]}")
    print(f"\nExplanation: {explanation}")

solve_medical_question()