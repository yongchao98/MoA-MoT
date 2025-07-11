def solve_disneyization_question():
    """
    Identifies the correct characteristics of Disneyfication based on Alan Bryman's theory.
    """
    # Alan Bryman's four main characteristics of Disneyization of Society (2004)
    bryman_characteristics = [
        "theming",
        "hybrid consumption",
        "merchandising",
        "performative labor"
    ]

    # The provided answer choices
    answer_choices = {
        "A": ["hybrid consumption", "merchandising"],
        "B": ["performative labor", "sanitization"],
        "C": ["trivialization", "theming"],
        "D": ["sanitization", "trivialization"],
        "E": ["disneyfication", "disneyization"],
        "F": ["mcdonaldization", "disneyization"],
        "G": ["theming", "performative labor"]
    }

    print("Alan Bryman's four key characteristics of Disneyization are:")
    for i, char in enumerate(bryman_characteristics, 1):
        print(f"{i}. {char.title()}")

    print("\nEvaluating the choices:")
    correct_choice = None
    for choice, items in answer_choices.items():
        # Check if both items in the choice are in Bryman's list of characteristics
        if all(item in bryman_characteristics for item in items):
            correct_choice = choice
            # While both A and G are technically correct, G pairs the foundational concept
            # of theming with the distinct human element of performative labor.
            # We will select G as the best representation.
            if choice == "G":
                break
    
    if correct_choice:
        print(f"\nChoice {correct_choice} contains two of Bryman's core characteristics:")
        chosen_characteristics = answer_choices[correct_choice]
        print(f"Characteristic 1: {chosen_characteristics[0]}")
        print(f"Characteristic 2: {chosen_characteristics[1]}")
    else:
        print("\nCould not find a choice with two of Bryman's primary characteristics.")

solve_disneyization_question()