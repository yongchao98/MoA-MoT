def identify_rocks():
    """
    Identifies the rock specimens based on their visual characteristics.
    """
    rock_identifications = {
        'H': 'Scoria - A reddish, vesicular (bubbly) extrusive igneous rock formed from gas-rich lava.',
        'I': 'Volcanic Breccia/Tuff - Composed of angular volcanic fragments and ash cemented together.',
        'J': 'Pumice - A light-colored, highly vesicular and lightweight volcanic rock.',
        'K': 'Pegmatite - A very coarse-grained intrusive igneous rock with large, interlocking crystals.'
    }

    print("Rock Specimen Identifications:")
    for specimen, description in rock_identifications.items():
        print(f"Specimen {specimen}: {description}")

identify_rocks()