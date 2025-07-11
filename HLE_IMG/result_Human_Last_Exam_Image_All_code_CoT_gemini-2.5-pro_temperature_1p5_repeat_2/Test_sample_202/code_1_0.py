def identify_rocks():
    """
    Identifies the rock specimens based on their visual characteristics.
    """
    specimens = {
        'H': 'Scoria - A vesicular (holey) extrusive igneous rock. The reddish color is due to iron oxidation.',
        'I': 'Volcanic Breccia - A rock composed of angular fragments of volcanic rock cemented together.',
        'J': 'Pumice - A very light-colored, extremely porous (vesicular) volcanic rock formed during explosive eruptions.',
        'K': 'Pegmatite - An intrusive igneous rock with very large, interlocking crystals, indicating slow cooling deep underground.'
    }

    print("Rock Specimen Identifications:")
    for specimen_id, description in specimens.items():
        print(f"Specimen {specimen_id}: {description}")

identify_rocks()