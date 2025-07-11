def identify_rocks():
    """
    Identifies the rock specimens shown in the image.
    """
    specimens = {
        'H': 'Scoria',
        'I': 'Concrete',
        'J': 'Pumice',
        'K': 'Gabbro'
    }

    print("Rock Specimen Identifications:")
    for letter, name in specimens.items():
        print(f"Specimen {letter} is {name}.")

identify_rocks()