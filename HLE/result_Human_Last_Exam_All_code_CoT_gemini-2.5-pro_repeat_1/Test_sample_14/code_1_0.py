def find_stunning_flyer():
    """
    Finds the T1 flying unit in BAR that can shoot and stun.
    """
    # A simplified database of T1 flying units in BAR.
    # Note: All units listed are from T1 air factories.
    t1_air_units = [
        {
            'name': 'Freedom Fighter',
            'faction': 'Armada',
            'type': 'Fighter',
            'abilities': ['shoot', 'anti-air']
        },
        {
            'name': 'Thunder',
            'faction': 'Armada',
            'type': 'Bomber',
            'abilities': ['shoot', 'bomb']
        },
        {
            'name': 'Vamp',
            'faction': 'Cortex',
            'type': 'Fighter',
            'abilities': ['shoot', 'anti-air']
        },
        {
            'name': 'Liche',
            'faction': 'Cortex',
            'type': 'Bomber',
            'abilities': ['shoot', 'bomb']
        },
        {
            'name': 'Instigator',
            'faction': 'Cortex',
            'type': 'Gunship',
            'abilities': ['shoot', 'stun', 'emp']
        }
    ]

    # Search for the unit that meets the criteria
    found_unit = None
    for unit in t1_air_units:
        # A unit is considered to be able to "shoot" if it has a weapon.
        # We check for the specific 'stun' ability.
        if 'shoot' in unit['abilities'] and 'stun' in unit['abilities']:
            found_unit = unit
            break

    if found_unit:
        print(f"The T1 flying unit that can shoot and stun enemies is the: {found_unit['name']}")
        print(f"It is a Cortex {found_unit['type']} that fires EMP rockets, which stun targets upon impact.")
    else:
        print("No T1 flying unit with the ability to stun was found.")

if __name__ == '__main__':
    find_stunning_flyer()