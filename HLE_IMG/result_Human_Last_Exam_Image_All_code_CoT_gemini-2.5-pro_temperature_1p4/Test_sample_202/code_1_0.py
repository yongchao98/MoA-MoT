def identify_rocks():
    """
    Identifies the rock specimens based on their visual characteristics.
    """
    specimens = {
        'H': 'Scoria',
        'I': 'Volcanic Breccia',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for label, rock_type in specimens.items():
        print(f"Specimen {label}: {rock_type}")

identify_rocks()