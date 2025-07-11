def identify_rocks():
    """
    This function identifies the rock specimens based on their visual characteristics.
    """
    identifications = {
        'H': 'Scoria',
        'I': 'Tuff',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for specimen, rock_type in identifications.items():
        print(f"Specimen {specimen}: {rock_type}")

identify_rocks()