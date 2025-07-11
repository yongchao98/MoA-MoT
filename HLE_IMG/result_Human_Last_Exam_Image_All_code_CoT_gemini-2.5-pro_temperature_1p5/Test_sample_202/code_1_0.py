def identify_rocks():
    """
    Identifies the rock specimens based on their visual characteristics.
    """
    specimens = {
        'H': 'Scoria - A reddish, vesicular (porous) volcanic rock.',
        'I': 'Volcanic Breccia / Tuff - A rock composed of angular volcanic fragments cemented together.',
        'J': 'Pumice - A light-colored, highly vesicular volcanic rock that is very lightweight.',
        'K': 'Pegmatite - A coarse-grained igneous rock with very large, interlocking crystals.'
    }

    print("Rock Specimen Identifications:")
    for label, description in specimens.items():
        print(f"Specimen {label}: {description}")

identify_rocks()