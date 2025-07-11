def identify_rocks():
    """
    Identifies the rock specimens based on their visual characteristics.
    """
    specimens = {
        'H': 'Scoria',
        'I': 'Breccia',
        'J': 'Pumice',
        'K': 'Granite/Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for letter, name in specimens.items():
        print(f"Specimen {letter}: {name}")

identify_rocks()