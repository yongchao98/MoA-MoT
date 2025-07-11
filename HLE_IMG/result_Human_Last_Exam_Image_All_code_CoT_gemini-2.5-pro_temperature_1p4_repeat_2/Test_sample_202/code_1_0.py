def identify_rocks():
    """
    Identifies the rock specimens based on visual characteristics.
    """
    specimens = {
        'H': 'Scoria. This is a vesicular (holey) extrusive igneous rock. The reddish color comes from the oxidation of iron as it cooled.',
        'I': 'Volcanic Tuff/Breccia. This rock is composed of consolidated volcanic ash and larger fragments (clasts) ejected during a volcanic eruption.',
        'J': 'Pumice. This is a light-colored, extremely porous (vesicular) volcanic rock formed during an explosive eruption. It has a low density and can often float on water.',
        'K': 'Pegmatite. This is a coarse-grained intrusive igneous rock. The large, visible crystals indicate it cooled very slowly deep within the Earth.'
    }

    print("Rock Specimen Identifications:")
    for specimen, description in specimens.items():
        print(f"Specimen {specimen}: {description}")

identify_rocks()