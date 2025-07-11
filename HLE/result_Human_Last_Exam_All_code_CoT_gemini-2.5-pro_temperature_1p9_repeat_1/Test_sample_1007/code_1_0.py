def find_ballet_step():
    """
    Analyzes classical ballet steps to find which one ends in the same position it started.
    """
    # Define the starting and ending positions for each step.
    # "Position A" can represent right foot front in fifth position.
    # "Position B" can represent left foot front in fifth position.
    # The key is to see if the end position is the same as the start, not what the position is.
    ballet_steps = [
        {
            "name": "A. Entrechat six",
            "start_position": "Position A",
            "end_position": "Position A",
            "explanation": "A jump with six crossings, landing in the original starting position."
        },
        {
            "name": "B. Échappé battu changé",
            "start_position": "Position A",
            "end_position": "Position B",
            "explanation": "The 'changé' (changed) indicates the feet switch places upon closing."
        },
        {
            "name": "C. Assemblé",
            "start_position": "Position A",
            "end_position": "Position B",
            "explanation": "Typically changes the feet (e.g., an assemblé dessus brings the back foot to the front)."
        },
        {
            "name": "D. Glissade derrière",
            "start_position": "Position A",
            "end_position": "Position B",
            "explanation": "A glide where the starting front foot finishes in the back."
        },
        {
            "name": "E. Gargouillade",
            "start_position": "Position A",
            "end_position": "Position B",
            "explanation": "A complex step that most often lands with the feet changed."
        }
    ]

    print("Analyzing ballet steps to find which ends in the same position it started from:\n")

    correct_answer = None
    for step in ballet_steps:
        is_same_position = step["start_position"] == step["end_position"]
        print(f"Step: {step['name']}")
        print(f"  - Starts in: {step['start_position']}")
        print(f"  - Ends in:   {step['end_position']}")
        print(f"  - Is the ending position the same as the starting position? {is_same_position}")
        print(f"  - Reason: {step['explanation']}")
        print("-" * 20)
        if is_same_position:
            correct_answer = step['name']

    if correct_answer:
        print(f"\nConclusion: The step that ends in the same position it started is {correct_answer}.")
    else:
        print("\nConclusion: No step was found that ends in the same position.")


find_ballet_step()