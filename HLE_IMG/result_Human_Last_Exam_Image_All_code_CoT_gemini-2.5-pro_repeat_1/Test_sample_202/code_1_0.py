def identify_rocks():
    """
    Identifies the rock specimens shown in the image.
    """
    specimens = {
        'H': 'Scoria',
        'I': 'Breccia',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for label, rock_type in specimens.items():
        print(f"Specimen {label}: {rock_type}")

identify_rocks()