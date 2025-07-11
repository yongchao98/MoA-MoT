def find_vaganova_arabesques():
    """
    Identifies Vaganova arabesques based on the relationship
    between the forward arm and the lifted leg.
    """
    # In Vaganova technique:
    # 1st Arabesque: Arm on the same side as the lifted leg is forward.
    # 2nd Arabesque: Arm on the opposite side as the lifted leg is forward.
    # 3rd Arabesque: Arm on the opposite side as the lifted leg is forward (arms are extended from 2nd position, one forward, one to the side, maintaining the line of the shoulders). Some descriptions vary, but the primary forward arm is opposite the lifted leg.
    # 4th Arabesque: Arm on the same side as the lifted leg is forward, with the body in a strong twist.
    
    arabesques_data = {
        'First':  {'number': 1, 'relation': 'same'},
        'Second': {'number': 2, 'relation': 'opposite'},
        'Third':  {'number': 3, 'relation': 'opposite'},
        'Fourth': {'number': 4, 'relation': 'same'}
    }

    matching_arabesques = []
    
    # Find arabesques where the forward arm is opposite the lifted leg
    for name, properties in arabesques_data.items():
        if properties['relation'] == 'opposite':
            matching_arabesques.append((name, properties['number']))

    print("The Vaganova arabesques where the forward arm is on the opposite side as the lifted leg are:")
    
    # Sort for consistent output and print each number
    matching_arabesques.sort(key=lambda x: x[1])
    
    final_numbers = []
    for name, number in matching_arabesques:
        print(f"- {name} Arabesque")
        final_numbers.append(str(number))

    # Per the instructions, showing each number in a "final equation"
    # which we interpret as listing the identified numbers.
    equation_str = " and ".join(final_numbers)
    print(f"\nThe corresponding arabesque numbers are {equation_str}.")
    print("\nThese two are the Second and Third arabesques, which corresponds to option B.")

find_vaganova_arabesques()