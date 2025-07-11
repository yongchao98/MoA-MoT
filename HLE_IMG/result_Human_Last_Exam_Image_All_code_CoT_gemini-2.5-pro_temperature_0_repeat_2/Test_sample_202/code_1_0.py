def identify_rocks():
    """
    Identifies the rock specimens shown in the image.
    """
    specimens = {
        'H': 'Scoria',
        'I': 'Tuff (or Volcanic Breccia)',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for letter, rock_type in specimens.items():
        print(f"Specimen {letter}: {rock_type}")

identify_rocks()