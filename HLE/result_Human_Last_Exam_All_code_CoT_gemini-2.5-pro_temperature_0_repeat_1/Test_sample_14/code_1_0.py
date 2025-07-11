def find_stunning_t1_flyer():
    """
    This function identifies flying units from T1 buildings in BAR
    that can attack and stun enemies.
    """
    # Data representing a selection of T1 flying units in BAR
    # Attributes: name, type, tier, can_attack, abilities
    t1_air_units = [
        {'name': 'Fink', 'type': 'air', 'tier': 1, 'can_attack': False, 'abilities': ['scout']},
        {'name': 'Vulture', 'type': 'air', 'tier': 1, 'can_attack': True, 'abilities': ['bomber']},
        {'name': 'Hawk', 'type': 'air', 'tier': 1, 'can_attack': True, 'abilities': ['stun']},
        {'name': 'Gnat', 'type': 'air', 'tier': 1, 'can_attack': True, 'abilities': []},
        {'name': 'Thunder', 'type': 'air', 'tier': 1, 'can_attack': True, 'abilities': ['bomber']},
    ]

    print("Searching for a T1 flying unit that can shoot and stun...")
    
    found_unit = None
    for unit in t1_air_units:
        # Check if the unit meets all criteria
        is_t1 = (unit['tier'] == 1)
        is_air = (unit['type'] == 'air')
        can_shoot = unit['can_attack']
        can_stun = ('stun' in unit['abilities'])

        if is_t1 and is_air and can_shoot and can_stun:
            found_unit = unit['name']
            break # Stop searching once we find a match

    if found_unit:
        print(f"Found a match: The {found_unit}.")
        print(f"This unit is from tier {1}, can fly, shoot, and has the 'stun' ability.")
    else:
        print("No T1 flying unit found with the ability to shoot and stun.")

# Execute the function to find and print the answer
find_stunning_t1_flyer()