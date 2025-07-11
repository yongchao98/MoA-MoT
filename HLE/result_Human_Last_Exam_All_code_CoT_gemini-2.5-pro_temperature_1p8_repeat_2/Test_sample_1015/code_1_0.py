def solve_ballet_question():
    """
    Analyzes and answers a multiple-choice question about ballet technique.
    """
    question = "What is the difference between a cambré derrière in the Vaganova and Balanchine methods?"
    
    options = {
        'A': 'Arm placement during allongé',
        'B': 'Placement of hip',
        'C': 'Speed',
        'D': 'Degree of backbend',
        'E': 'Placement of head'
    }

    correct_answer_key = 'B'
    
    explanation = (
        "In a Vaganova cambré derrière, the technique emphasizes keeping the hips square and level, "
        "directly over the supporting foot. The bend originates from the upper spine, creating a clean, "
        "supported arc. The body acts as a single, controlled unit.\n\n"
        "In the Balanchine method, a key stylistic difference is the treatment of the hip. The hip of the "
        "supporting leg is often pushed diagonally forward and out to the side. This 'released' hip "
        "creates a more extreme, dynamic line and a visually distinct C-curve shape. It's a foundational "
        "difference in mechanics that alters the entire look of the step.\n\n"
        "While other elements like head or arm placement can also vary, the fundamental "
        "technical distinction that defines the style difference is the placement of the hip."
    )

    print(f"Question: {question}\n")
    print("Choices:")
    for key, value in options.items():
        print(f"  {key}. {value}")

    print("\n--------------------------")
    print(f"ANALYSIS AND ANSWER")
    print("--------------------------")
    print(f"The correct option is: {correct_answer_key}")
    print(f"The correct answer is: '{options[correct_answer_key]}'")
    print("\nExplanation:")
    print(explanation)

solve_ballet_question()