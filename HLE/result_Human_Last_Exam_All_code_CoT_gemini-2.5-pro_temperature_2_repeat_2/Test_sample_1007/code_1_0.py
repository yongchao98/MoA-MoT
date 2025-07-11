def find_ballet_step():
    """
    Analyzes classical ballet steps to find which one ends in the same
    leg position as it started.
    """

    # A dictionary representing the ballet steps and whether they change leg position.
    # True: Starting and ending leg positions are different.
    # False: Starting and ending leg positions are the same.
    ballet_steps = {
        "A. Entrechat six": {
            "changes_position": True,
            "description": "A jump with six crossings, landing with the opposite foot in front."
        },
        "B. Échappé battu changé": {
            "changes_position": True,
            "description": "A jump opening from fifth to second and closing back to fifth with the feet changed."
        },
        "C. Assemblé": {
            "changes_position": False,
            "description": "A jump where the legs are assembled in the air to land in the original starting position."
        },
        "D. Glissade derrière": {
            "changes_position": True,
            "description": "A gliding step that typically results in the feet changing their front/back position."
        },
        "E. Gargouillade": {
            "changes_position": True,
            "description": "A complex jump where the dancer lands with the feet switched from the starting position."
        }
    }

    print("Analyzing ballet steps to find the one with the same start and end position...")
    correct_answer = None
    for step, properties in ballet_steps.items():
        if not properties["changes_position"]:
            correct_answer = step
            break # Found the answer, no need to check further.

    if correct_answer:
        print(f"\nThe correct answer is: {correct_answer}")
        print(f"Reason: {ballet_steps[correct_answer]['description']}")
    else:
        print("No step found that meets the criteria.")

find_ballet_step()