def identify_rocks():
    """
    Identifies the rock specimens based on visual characteristics.
    """
    specimens = {
        'H': 'Scoria',
        'I': 'Tuff',
        'J': 'Pumice',
        'K': 'Pegmatite'
    }

    print("Rock Specimen Identifications:")
    for specimen_id, rock_name in specimens.items():
        print(f"Specimen {specimen_id} is {rock_name}.")

identify_rocks()