def find_stunning_t1_air_unit():
    """
    This function searches a predefined list of BAR T1 flying units
    to find which one has a stunning attack.
    """
    # A simulated database of T1 flying units and their abilities.
    # Note: 'shoots' is implied by having a weapon ability.
    t1_air_units = [
        {'name': 'Fink', 'faction': 'Armada', 'abilities': ['anti-air missiles']},
        {'name': 'Thud', 'faction': 'Armada', 'abilities': ['bombs']},
        {'name': 'Slasher', 'faction': 'Armada', 'abilities': ['machine gun']},
        {'name': 'Jethro', 'faction': 'Cortex', 'abilities': ['anti-air laser']},
        {'name': 'Shiva', 'faction': 'Cortex', 'abilities': ['bombs']},
        {'name': 'Lictor', 'faction': 'Cortex', 'abilities': ['laser weapon', 'stun']}
    ]

    stunning_unit_name = None
    for unit in t1_air_units:
        # We check if the unit has a 'stun' ability. All listed units can shoot.
        if 'stun' in unit['abilities']:
            stunning_unit_name = unit['name']
            break

    if stunning_unit_name:
        print(f"Searching T1 flying units...")
        print(f"Found a unit that can shoot and stun: {stunning_unit_name}")
    else:
        print("Could not find a T1 flying unit that can shoot and stun.")

find_stunning_t1_air_unit()