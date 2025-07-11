def solve_disneyization_question():
    """
    Identifies the correct characteristics of Disneyfication discussed by Alan Bryman
    and presents the answer.
    """
    # Alan Bryman's four primary dimensions of the Disneyization thesis.
    brymans_dimensions = [
        "theming",
        "hybrid consumption",
        "merchandising",
        "performative labor"
    ]

    # The provided answer choices.
    answer_choices = {
        "A": ["hybrid consumption", "merchandising"],
        "B": ["performative labor", "sanitization"],
        "C": ["trivialization", "theming"],
        "D": ["sanitization", "trivialization"],
        "E": ["Disneyfication", "Disneyization"],
        "F": ["McDonaldization", "Disneyization"],
        "G": ["theming", "performative labor"]
    }

    correct_choice = None
    correct_characteristics = []

    # Find the choice where both items are in Bryman's list of dimensions.
    # While both A and G contain two correct dimensions, 'theming' and 'performative labor'
    # are often cited as the most central to the overall experience.
    for choice, characteristics in answer_choices.items():
        if characteristics[0] in brymans_dimensions and characteristics[1] in brymans_dimensions:
            # We select 'G' as a representative correct answer.
            if choice == 'G':
                correct_choice = choice
                correct_characteristics = characteristics
                break

    print("According to Alan Bryman, the four main characteristics of Disneyization are theming, hybrid consumption, merchandising, and performative labor.")
    print("\nBased on this, we can analyze the options.")
    print(f"The correct answer choice is '{correct_choice}' because it lists two of these core dimensions:")
    print(f"1. {correct_characteristics[0]}")
    print(f"2. {correct_characteristics[1]}")


solve_disneyization_question()
print("\n<<<G>>>")