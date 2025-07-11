def identify_rocks():
    """
    Identifies the rock specimens based on their visual characteristics.
    """
    # A dictionary to store the identification of each specimen.
    specimen_identifications = {
        'H': 'Scoria - A vesicular (holey) volcanic rock formed from gas-rich lava. The reddish color is due to iron oxidation.',
        'I': 'Breccia - A sedimentary or fault-related rock composed of broken, angular fragments of minerals or rock cemented together.',
        'J': 'Pumice - A very light and porous volcanic rock formed when gas-rich, frothy lava cools rapidly.',
        'K': 'Pegmatite - An intrusive igneous rock known for its very large, interlocking crystals that form from the slow cooling of magma.'
    }

    print("Rock Specimen Identifications:")
    print("=" * 30)
    for specimen, identification in specimen_identifications.items():
        print(f"Specimen {specimen}: {identification}")
    print("=" * 30)

# Run the identification function
identify_rocks()