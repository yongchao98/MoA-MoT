def solve_ballet_question():
    """
    Analyzes classical ballet steps to find which one starts and ends
    in the same leg position.
    """
    steps = [
        {
            "option": "A",
            "name": "Entrechat six",
            "description": "Starts from fifth position, beats the legs six times in the air, and lands in the same fifth position.",
            "ends_in_start_position": True
        },
        {
            "option": "B",
            "name": "Échappé battu changé",
            "description": "Starts in fifth position and lands in fifth position with the feet changed ('changé').",
            "ends_in_start_position": False
        },
        {
            "option": "C",
            "name": "Assemblé",
            "description": "A step where the working foot brushes and joins the other, typically changing the front foot upon landing in fifth.",
            "ends_in_start_position": False
        },
        {
            "option": "D",
            "name": "Glissade derrière",
            "description": "A gliding step where the front foot closes behind, changing the starting fifth position.",
            "ends_in_start_position": False
        },
        {
            "option": "E",
            "name": "Gargouillade",
            "description": "A complex jump that typically involves changing the feet upon landing.",
            "ends_in_start_position": False
        }
    ]

    correct_answer = None
    print("Analyzing each ballet step:")
    print("-" * 30)

    for step in steps:
        print(f"Option {step['option']}: {step['name']}")
        print(f"Description: {step['description']}")
        if step['ends_in_start_position']:
            correct_answer = step
            print("Result: Starts and ends in the SAME position.\n")
        else:
            print("Result: Starts and ends in a DIFFERENT position.\n")
    
    if correct_answer:
        print("-" * 30)
        print("Conclusion:")
        print(f"The correct answer is '{correct_answer['name']}', as it is the only step listed that concludes in the same position from which it began.")
        print(f"Final Answer Option: {correct_answer['option']}")

solve_ballet_question()