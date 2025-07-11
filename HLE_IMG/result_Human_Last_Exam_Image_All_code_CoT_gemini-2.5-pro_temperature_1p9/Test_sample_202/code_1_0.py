def identify_rocks():
    """
    Identifies the rock specimens shown in the image.
    """
    specimens = {
        'H': 'Scoria - A vesicular (holey) extrusive igneous rock.',
        'I': 'Breccia - A sedimentary rock composed of large, angular fragments cemented together.',
        'J': 'Pumice - A very light and porous volcanic rock that forms during explosive eruptions.',
        'K': 'Granite - A coarse-grained, intrusive igneous rock composed mainly of quartz and feldspar with minor micas and amphiboles.'
    }

    print("Rock Specimen Identifications:")
    for specimen, description in specimens.items():
        print(f"Specimen {specimen}: {description}")

identify_rocks()