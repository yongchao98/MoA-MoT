def find_stunning_flyer():
    """
    This function searches through a predefined list of BAR T1 flying units
    to find one that can both shoot and stun.
    """
    t1_flying_units = [
        {
            'name': 'Swift',
            'faction': 'Armada',
            'role': 'Fighter',
            'can_shoot': True,
            'can_stun': False
        },
        {
            'name': 'Hawk',
            'faction': 'Armada',
            'role': 'Gunship',
            'can_shoot': True,
            'can_stun': False
        },
        {
            'name': 'Thud',
            'faction': 'Cortex',
            'role': 'Gunship',
            'can_shoot': True,
            'can_stun': True
        },
        {
            'name': 'Vampire',
            'faction': 'Cortex',
            'role': 'Fighter',
            'can_shoot': True,
            'can_stun': False
        }
    ]

    for unit in t1_flying_units:
        # Check for a T1 unit that can both shoot and stun
        if unit['can_shoot'] and unit['can_stun']:
            print(f"The T1 flying unit that can shoot and stun is: {unit['name']}")
            return

find_stunning_flyer()