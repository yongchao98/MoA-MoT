def solve_arabesque_question():
    """
    Identifies Vaganova arabesques based on arm and leg position.
    In the Vaganova technique, each of the four basic arabesques has a specific
    rule for the placement of the arms relative to the legs. This function
    encodes those rules and finds the ones that match the specified criteria.
    """

    # Data describing the arm positions for the four Vaganova arabesques.
    # The 'arm_position' key describes the forward arm's relation to the lifted leg.
    arabesque_positions = [
        {'name': 'First', 'arm_position': 'opposite'},
        {'name': 'Second', 'arm_position': 'same'},
        {'name': 'Third', 'arm_position': 'opposite'},
        {'name': 'Fourth', 'arm_position': 'same'}
    ]

    # The condition from the question
    condition = 'opposite'

    # Find the arabesques that match the condition
    matching_positions = []
    for pos in arabesque_positions:
        if pos['arm_position'] == condition:
            matching_positions.append(pos['name'])

    print(f"The question asks for the arabesques where the front arm is on the '{condition}' side as the lifted leg.")
    print("Based on the rules of the Vaganova technique, the matching positions are:")
    print(" and ".join(matching_positions))


solve_arabesque_question()