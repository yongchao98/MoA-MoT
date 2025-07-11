def find_arabesques():
    """
    Identifies the Vaganova arabesques where the forward arm is on the
    opposite side of the lifted leg.
    """
    # In Vaganova terminology:
    # 'opposite' means the forward arm is on the same side as the *supporting* leg.
    # 'same' means the forward arm is on the same side as the *lifted* leg.
    arabesque_definitions = {
        "First": "opposite",
        "Second": "same",
        "Third": "opposite",
        "Fourth": "same"
    }

    matching_arabesques = []
    for position, arm_leg_relation in arabesque_definitions.items():
        if arm_leg_relation == "opposite":
            matching_arabesques.append(position)

    print("In the Vaganova technique, the arabesques where the dancer's forward arm is on the opposite side of the lifted leg are:")
    # The 'end' parameter prevents a newline, and the 'sep' parameter sets the separator.
    print(*matching_arabesques, sep=" and ")

find_arabesques()