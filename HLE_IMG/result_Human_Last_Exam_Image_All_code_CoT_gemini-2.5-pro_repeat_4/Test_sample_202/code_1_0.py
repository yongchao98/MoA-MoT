def identify_rocks():
    """
    Identifies the rock specimens based on visual characteristics.
    """
    specimens = {
        'H': 'Scoria',
        'I': 'Breccia',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for letter, rock_type in specimens.items():
        print(f"Specimen {letter} is {rock_type}.")

identify_rocks()