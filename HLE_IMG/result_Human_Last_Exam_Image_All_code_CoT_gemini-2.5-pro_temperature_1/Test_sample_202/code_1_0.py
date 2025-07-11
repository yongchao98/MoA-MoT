def identify_rocks():
    """
    Identifies the rock specimens based on their visual characteristics.
    """
    specimens = {
        'H': 'Scoria - A vesicular (holey) volcanic rock, typically dark in color (brown, black, or red).',
        'I': 'Breccia - A rock composed of broken, angular fragments of minerals or rock cemented together by a fine-grained matrix.',
        'J': 'Pumice - A very light and porous volcanic rock formed when a gas-rich froth of glassy lava solidifies rapidly.',
        'K': 'Pegmatite - A coarse-grained intrusive igneous rock with very large, interlocking crystals, typically of a composition similar to granite.'
    }

    print("Rock Specimen Identifications:")
    print("----------------------------")
    for letter, description in specimens.items():
        print(f"Specimen {letter}: {description}")

identify_rocks()