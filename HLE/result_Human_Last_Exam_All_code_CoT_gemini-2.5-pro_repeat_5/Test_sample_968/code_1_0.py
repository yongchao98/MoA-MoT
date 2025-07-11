def find_matching_arabesques():
    """
    Identifies Vaganova arabesques where the forward arm is on the opposite side of the lifted leg.
    """

    # In ballet, 'opposition' means the arm and leg are on opposite sides.
    # We will define the four Vaganova arabesques based on this principle.
    # A value of 'True' means the forward arm is opposite the lifted leg.
    # A value of 'False' means it's on the same side.
    vaganova_arabesques = {
        1: {'is_opposition': True, 'name': 'First Arabesque'},
        2: {'is_opposition': False, 'name': 'Second Arabesque'},
        3: {'is_opposition': True, 'name': 'Third Arabesque'},
        4: {'is_opposition': False, 'name': 'Fourth Arabesque'}
    }

    matching_numbers = []

    print("Analyzing the Vaganova Arabesques:")
    print("-----------------------------------")
    for number, details in vaganova_arabesques.items():
        if details['is_opposition']:
            print(f"- {details['name']} ({number}): The forward arm is OPPOSITE the lifted leg. This matches.")
            matching_numbers.append(number)
        else:
            print(f"- {details['name']} ({number}): The forward arm is on the SAME SIDE as the lifted leg. This does not match.")
    
    print("\n-----------------------------------")
    print("Final Result:")
    print("The two types of arabesque that fit the description are the First and Third.")
    
    # As requested, printing the numbers for the final equation.
    print(f"The final equation consists of the numbers for the matching arabesques:")
    print(f"Number 1: {matching_numbers[0]}")
    print(f"Number 2: {matching_numbers[1]}")

find_matching_arabesques()