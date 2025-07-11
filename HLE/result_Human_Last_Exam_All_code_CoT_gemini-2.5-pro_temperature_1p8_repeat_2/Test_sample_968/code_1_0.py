def find_arabesques():
    """
    Identifies the Vaganova arabesques where the forward arm is
    on the opposite side of the lifted leg.
    """
    # In Vaganova technique:
    # 'opposite' means forward arm is opposite to the lifted leg (e.g., right leg up, left arm forward).
    # 'same' means forward arm is on the same side as the lifted leg (e.g., right leg up, right arm forward).
    vaganova_rules = {
        'First': 'opposite',
        'Second': 'same',
        'Third': 'opposite',
        'Fourth': 'same'
    }

    print("Analyzing the Vaganova arabesque positions:")
    print(f" - First Arabesque: The forward arm is on the '{vaganova_rules['First']}' side as the lifted leg.")
    print(f" - Second Arabesque: The forward arm is on the '{vaganova_rules['Second']}' side as the lifted leg.")
    print(f" - Third Arabesque: The forward arm is on the '{vaganova_rules['Third']}' side as the lifted leg.")
    print(f" - Fourth Arabesque: The forward arm is on the '{vaganova_rules['Fourth']}' side as the lifted leg.")

    matching_arabesques = []
    for arabesque, arm_position in vaganova_rules.items():
        if arm_position == 'opposite':
            matching_arabesques.append(arabesque)
            
    print("\nThe two arabesques with the arm extended on the opposite side of the lifted leg are:")
    # Per the instructions, printing the "numbers" in the final equation.
    print(f"{matching_arabesques[0]} and {matching_arabesques[1]}")

find_arabesques()