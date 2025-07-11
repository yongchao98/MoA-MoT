def find_photochemical_synthesizers():
    """
    Identifies taxa that perform photochemical synthesis from a predefined list.

    The function evaluates a list of species to determine if they use light
    energy to synthesize compounds (e.g., photosynthesis, bacteriorhodopsin-based
    phototrophy) as a normal metabolic process, ignoring symbiotic relationships.
    """

    # A dictionary mapping index to a tuple: (name, performs_photochemical_synthesis)
    # The boolean is determined based on the organism's known metabolic pathways.
    taxa_data = {
        1: ("Acanthella cavernosa", False),           # Animal (Sponge)
        2: ("Gloeochaete wittrockiana", True),      # Alga (Photosynthesis)
        3: ("Homo sapiens", False),                 # Animal (Human)
        4: ("Riftia pachyptila", False),            # Animal (Tube Worm)
        5: ("Halapricum salinum", True),            # Archaea (Phototrophy)
        6: ("Aphanothece castagnei", True),         # Cyanobacteria (Photosynthesis)
        7: ("Baileya pleniradiata", True),          # Plant (Photosynthesis)
        8: ("Acanthella pulchra", False),           # Animal (Sponge)
        9: ("Ectothiorhodosinus mongolicus", True), # Bacteria (Anoxygenic Photosynthesis)
        10: ("Chlorobaculum tepidum", True),        # Bacteria (Anoxygenic Photosynthesis)
        11: ("Stygichthys typhlops", False),        # Animal (Fish)
        12: ("Gemmatimonas phototrophica", True),   # Bacteria (Anoxygenic Photosynthesis)
        13: ("Myonera garretti", False)              # Animal (Bivalve)
    }

    phototrophic_indices = []
    for index, (name, is_phototrophic) in taxa_data.items():
        if is_phototrophic:
            phototrophic_indices.append(str(index))

    if not phototrophic_indices:
        print("none")
    else:
        # Join the collected indices with commas
        result_string = ",".join(phototrophic_indices)
        print(result_string)

# Execute the function to find and print the result
find_photochemical_synthesizers()