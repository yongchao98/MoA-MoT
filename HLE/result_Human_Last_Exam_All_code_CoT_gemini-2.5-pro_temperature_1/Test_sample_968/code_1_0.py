def find_matching_arabesques():
    """
    Identifies Vaganova arabesques where the forward arm is on the
    opposite side of the lifted leg.
    """
    # Define the Vaganova arabesques.
    # The value represents the position of the forward arm relative to the lifted leg.
    vaganova_arabesques = {
        "First": "opposite",
        "Second": "same",
        "Third": "opposite",
        "Fourth": "same"
    }

    matching_types = []
    for arabesque_type, arm_position in vaganova_arabesques.items():
        if arm_position == "opposite":
            matching_types.append(arabesque_type)

    if len(matching_types) == 2:
        print(f"The two types of arabesque that fit the description are the {matching_types[0]} and {matching_types[1]}.")
        print("This corresponds to answer choice C.")
    else:
        print("Could not determine the two matching arabesque types.")

find_matching_arabesques()