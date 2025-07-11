import textwrap

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

    explanation = """
    The most significant and defining technical difference between a cambré derrière in the Vaganova and Balanchine methods is the placement of the hip.

    - **Vaganova Method:** Emphasizes a strict, classical form where the hips and shoulders are kept square to the front. The backbend is initiated from the upper back, maintaining a stable and aligned torso. The integrity of this 'box' alignment is paramount.

    - **Balanchine Method:** Often encourages a more dynamic and extreme aesthetic. To achieve a deeper and more sweeping cambré derrière, the Balanchine technique involves releasing the hip of the working leg, allowing it to press forward. This breaks the square alignment of the torso but creates a visually dramatic line.

    While other factors like the degree of the backbend (D) and the placement of the head (E) also differ, they are direct consequences of this fundamental difference in hip placement. Releasing the hip (B) is what allows for the deeper backbend. Therefore, the placement of the hip is the root distinction.
    """

    print(question)
    for key, value in options.items():
        print(f"{key}. {value}")
    
    print("\n--- Explanation ---")
    print(textwrap.dedent(explanation).strip())
    print("-------------------")
    print(f"\nThe correct answer is B: {options[correct_answer_key]}.")


solve_ballet_question()
<<<B>>>