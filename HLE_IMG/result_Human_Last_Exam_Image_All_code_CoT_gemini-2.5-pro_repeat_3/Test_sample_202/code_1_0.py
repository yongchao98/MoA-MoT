def identify_rocks():
    """
    This function identifies the rock specimens shown in the image.
    """
    specimens = {
        'H': 'Scoria',
        'I': 'Tuff',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for letter, rock_name in specimens.items():
        print(f"Specimen {letter} is {rock_name}.")

identify_rocks()