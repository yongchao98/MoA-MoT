def find_ballet_step():
    """
    Analyzes classical ballet steps to identify which one has the same
    starting and ending leg position.
    """
    # A database of ballet steps and their properties.
    # The key 'ends_same_position' is True if the starting and ending
    # positions of the feet are identical (e.g., right foot front in fifth
    # at the start and right foot front in fifth at the end).
    ballet_steps = {
        "A. Entrechat six": {
            "description": "Starts in fifth, jumps, beats six times, and lands in fifth with the feet changed.",
            "ends_same_position": False
        },
        "B. Échappé battu changé": {
            "description": "Starts in fifth, jumps to second, beats, and lands back in fifth with the feet changed.",
            "ends_same_position": False
        },
        "C. Assemblé": {
            "description": "Starts from fifth, one leg extends, pushes off the supporting leg, and 'assembles' the legs back into the original starting fifth position upon landing.",
            "ends_same_position": True
        },
        "D. Glissade derrière": {
            "description": "A traveling step starting from fifth. The back foot glides out and the front foot closes behind, changing the feet.",
            "ends_same_position": False
        },
        "E. Gargouillade": {
            "description": "A decorative step involving circular movements of the legs in the air, landing with the feet changed.",
            "ends_same_position": False
        }
    }

    correct_answer_key = None
    explanation = ""

    # Iterate through the steps to find the correct one
    for step, properties in ballet_steps.items():
        if properties["ends_same_position"]:
            correct_answer_key = step
            explanation = properties["description"]
            break

    if correct_answer_key:
        print(f"The ballet step with the same starting and ending leg position is:")
        print(f"{correct_answer_key}\n")
        print("Explanation:")
        print(explanation)
    else:
        print("Could not determine the correct answer based on the provided data.")

# Execute the function to find and print the answer.
find_ballet_step()
<<<C>>>