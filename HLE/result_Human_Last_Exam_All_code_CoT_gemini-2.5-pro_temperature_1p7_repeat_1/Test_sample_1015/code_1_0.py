def explain_ballet_technique_difference():
    """
    Explains the difference between a Vaganova and Balanchine cambré derrière.
    """
    question = "What is the difference between a cambré derrière in the Vaganova and Balanchine methods?"
    
    choices = {
        'A': 'Arm placement during allongé',
        'B': 'Placement of hip',
        'C': 'Speed',
        'D': 'Degree of backbend',
        'E': 'Placement of head'
    }

    explanation = """
The most defining technical difference between a cambré derrière in the Vaganova and Balanchine methods is the placement and action of the hips.

1.  **Vaganova Method:** In the Russian Vaganova method, the dancer keeps the hips square and stable over the supporting legs. The bend originates from the upper back (thoracic spine), with a feeling of lifting the torso up and back. The hips do not press forward; this ensures stability and control.

2.  **Balanchine Method:** The American Balanchine style is characterized by a dynamic push of the hips forward as the upper body bends backward. This creates a visually deeper, more extreme arc and contributes to the speed and 'off-balance' aesthetic for which Balanchine is known.

While other elements like head placement (E) and the degree of the backbend (D) also differ, they are largely a result of the fundamental difference in hip placement (B). Therefore, the placement of the hip is the core distinction.
"""

    correct_answer_letter = 'B'

    print(question)
    print("-" * 20)
    print("Answer Choices:")
    for letter, description in choices.items():
        print(f"{letter}. {description}")
    print("-" * 20)
    print("Explanation:")
    print(explanation)
    print(f"The most accurate answer is B: {choices[correct_answer_letter]}")

explain_ballet_technique_difference()