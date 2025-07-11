def find_matching_arabesques():
    """
    Identifies Vaganova arabesques based on arm and leg position.

    The function checks which arabesques have the forward-extended arm
    on the opposite side of the body as the lifted leg.
    """

    # Define the relationship between the forward arm and the lifted leg for each arabesque.
    # 'same': Forward arm is on the same side as the lifted leg.
    # 'opposite': Forward arm is on the opposite side of the lifted leg.
    vaganova_arabesques = {
        'First': 'opposite',
        'Second': 'same',
        'Third': 'opposite',
        'Fourth': 'same'
    }

    matching_types = []
    print("Analyzing Vaganova Arabesques:")
    for arabesque, arm_leg_relation in vaganova_arabesques.items():
        if arm_leg_relation == 'opposite':
            matching_types.append(arabesque)
            print(f"- {arabesque} Arabesque: The forward arm is on the OPPOSITE side of the lifted leg. (Match)")
        else:
            print(f"- {arabesque} Arabesque: The forward arm is on the SAME side of the lifted leg. (No Match)")

    print("\nConclusion:")
    print(f"The two types of arabesque that fit the criteria are the {matching_types[0]} and {matching_types[1]}.")
    print("This corresponds to answer choice C.")

find_matching_arabesques()