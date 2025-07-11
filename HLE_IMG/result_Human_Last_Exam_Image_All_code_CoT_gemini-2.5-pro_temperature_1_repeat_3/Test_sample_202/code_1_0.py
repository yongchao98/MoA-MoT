def identify_rocks():
    """
    Identifies the rock specimens based on visual characteristics.
    """
    specimens = {
        'H': 'Scoria',
        'I': 'Conglomerate',
        'J': 'Pumice',
        'K': 'Granite'
    }

    print("Rock Specimen Identifications:")
    for letter, rock_type in specimens.items():
        print(f"Specimen {letter}: {rock_type}")

identify_rocks()