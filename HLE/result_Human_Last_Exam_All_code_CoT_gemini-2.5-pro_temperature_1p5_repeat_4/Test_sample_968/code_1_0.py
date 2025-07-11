def find_matching_arabesques():
    """
    Identifies Vaganova arabesques based on arm and leg positions.
    """
    # Define the properties of the four Vaganova arabesques.
    # The key is the arabesque number, and the value describes the forward arm's
    # position relative to the lifted leg.
    # "opposite" means the arm on the same side as the *supporting* leg is forward.
    # "same" means the arm on the same side as the *lifted* leg is forward.
    vaganova_arabesques = {
        1: "opposite",
        2: "same",
        3: "opposite",
        4: "same"
    }

    matching_positions = []
    print("Analyzing Vaganova Arabesque Positions:")
    print("-----------------------------------------")
    for position, arm_leg_relation in vaganova_arabesques.items():
        if arm_leg_relation == "opposite":
            matching_positions.append(position)
            print(f"Position {position}: The forward arm is on the OPPOSITE side of the lifted leg. (Match)")
        else:
            print(f"Position {position}: The forward arm is on the SAME side as the lifted leg. (No Match)")
    
    print("\nConclusion:")
    print("The two types of arabesque that fit the description are the")
    # This line fulfills the requirement to output each number in the final equation/statement.
    print(f"first ({matching_positions[0]}) and third ({matching_positions[1]}).")


find_matching_arabesques()