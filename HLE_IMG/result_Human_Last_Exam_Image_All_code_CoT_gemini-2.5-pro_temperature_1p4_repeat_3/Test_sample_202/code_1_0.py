def identify_rocks():
    """
    Identifies the rock specimens based on their visual characteristics.
    """
    rock_identifications = {
        'H': 'Scoria',
        'I': 'Tuff',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for specimen, rock_type in rock_identifications.items():
        print(f"Specimen {specimen} is identified as: {rock_type}")

identify_rocks()