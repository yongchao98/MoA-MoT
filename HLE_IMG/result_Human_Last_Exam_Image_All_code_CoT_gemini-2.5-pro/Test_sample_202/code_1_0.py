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
    for letter, rock_type in specimens.items():
        print(f"Specimen {letter} is identified as: {rock_type}")

identify_rocks()