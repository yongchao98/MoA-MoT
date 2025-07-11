def identify_specimens():
    """
    Identifies the rock specimens based on their visual characteristics.
    """
    specimen_identifications = {
        'H': 'Scoria',
        'I': 'Tuff',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for specimen, rock_type in specimen_identifications.items():
        print(f"Specimen {specimen} is {rock_type}.")

identify_specimens()