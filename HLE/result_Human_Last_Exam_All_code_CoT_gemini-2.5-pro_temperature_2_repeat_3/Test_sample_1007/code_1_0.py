def find_ballet_step():
    """
    Analyzes classical ballet steps to find which one ends in the same
    leg position as it started.
    """
    # A dictionary to represent each ballet step and its properties.
    # 'changes_position' is True if the feet switch places (e.g., right-front to left-front).
    # 'changes_position' is False if the feet return to their original configuration.
    ballet_steps = [
        {
            "option": "A",
            "name": "Entrechat six",
            "description": "A jump with six leg crossings, landing in the original starting position.",
            "changes_position": False
        },
        {
            "option": "B",
            "name": "Échappé battu changé",
            "description": "A jump to second and back to fifth, with the feet changing places ('changé').",
            "changes_position": True
        },
        {
            "option": "C",
            "name": "Assemblé",
            "description": "A jump where the legs are assembled in the air, typically changing the foot position upon landing.",
            "changes_position": True
        },
        {
            "option": "D",
            "name": "Glissade derrière",
            "description": "A gliding step where the front foot closes behind, changing the foot position.",
            "changes_position": True
        },
        {
            "option": "E",
            "name": "Gargouillade",
            "description": "A complex jump involving ronds de jambe, not typically ending in the starting position.",
            "changes_position": True
        }
    ]

    print("Analyzing which ballet step starts and ends in the same leg position:")
    print("-" * 60)

    correct_step = None
    for step in ballet_steps:
        if not step["changes_position"]:
            correct_step = step
            break # Found the correct step

    if correct_step:
        print(f"The correct answer is '{correct_step['option']}. {correct_step['name']}'.")
        print(f"Reason: {correct_step['description']}")
    else:
        print("No step found that meets the criteria.")

find_ballet_step()