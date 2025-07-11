def find_ballet_step():
    """
    Analyzes classical ballet steps to find which one ends in the same
    leg position as it starts.
    """
    steps = {
        "A": {
            "name": "Entrechat six",
            "start_position": "Fifth, right foot front",
            "end_position": "Fifth, right foot front",
            "description": "A jump with three beats in the air, landing in the same fifth position."
        },
        "B": {
            "name": "Échappé battu changé",
            "start_position": "Fifth, right foot front",
            "end_position": "Fifth, left foot front",
            "description": "A jump from a closed to an open position and back to a closed position, changing which foot is in front."
        },
        "C": {
            "name": "Assemblé (dessus)",
            "start_position": "Fifth, left foot front",
            "end_position": "Fifth, right foot front",
            "description": "A jump where the working leg brushes out and 'assembles' in the air to land in front of the other."
        },
        "D": {
            "name": "Glissade derrière",
            "start_position": "Fifth, right foot front",
            "end_position": "Fifth, left foot front",
            "description": "A gliding step where the front foot ends up behind the other."
        },
        "E": {
            "name": "Gargouillade (dessus)",
            "start_position": "Fifth, right foot front",
            "end_position": "Fifth, left foot front",
            "description": "A complex jump similar to a pas de chat that results in changing the front foot."
        }
    }

    print("Analyzing which ballet step starts and ends in the same position...")
    print("-" * 60)

    correct_answer = None
    for key, step in steps.items():
        print(f"Step: {step['name']}")
        print(f"  Start: {step['start_position']}")
        print(f"  End:   {step['end_position']}")
        
        is_same = step["start_position"] == step["end_position"]
        print(f"  Positions are the same: {is_same}\n")
        
        if is_same:
            correct_answer = step['name']
            correct_key = key

    print("-" * 60)
    if correct_answer:
        print(f"Conclusion: {correct_answer} is the step that ends in the same position it starts from.")
        print(f"The correct answer choice is {correct_key}.")
    else:
        print("Could not find a step that meets the criteria in the provided list.")

if __name__ == "__main__":
    find_ballet_step()