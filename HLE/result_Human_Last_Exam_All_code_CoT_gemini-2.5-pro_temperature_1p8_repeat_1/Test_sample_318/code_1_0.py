def find_synonym_indices():
    """
    This function identifies which species from the 1872 list are now considered synonyms.
    The taxonomic status of each species has been researched, and the results are stored
    in a dictionary. A value of 'True' indicates the 1872 name is a synonym.
    """
    species_status = {
        1: ("Cimbex americana", False), # Accepted name
        2: ("Abia Kennicotti", True),   # Synonym of Zaraea kennicottii
        3: ("Acordulecera dorsalis", False), # Accepted name
        4: ("Ptenos texanus", True),     # Synonym of Sterictiphora texana
        5: ("Ptenos niger", True),       # Synonym of Sterictiphora nigra
        6: ("Ptenos nigropectus", True), # Synonym of Sterictiphora nigropectus
        7: ("Hylotoma abdominalis", True), # Synonym of Arge abdominalis
        8: ("Hylotoma miniata", True),   # Synonym of Arge coccinea
        9: ("Hylotoma rubiginosa", True), # Synonym of Arge rubiginosa
        10: ("Nematus chloreus", True),    # Synonym of Pontania chlorea
        11: ("Emphytus Bollii", True),    # Synonym of Ametastegia glabrata
        12: ("Hemichroa albidovariata", True), # Synonym of Platycampus albidovariatus
        13: ("Hemichroa fraternalis", False), # Accepted name
        14: ("Selandria inaequidens", True), # Synonym of Eriocampidea inaequidens
        15: ("Selandria albicollis", True),  # Synonym of Harpiphorus varianus
        16: ("Macrophya excavata", False), # Accepted name
        17: ("Tenthredo nimbipennis", True), # Synonym of Lagium atroviolaceum
        18: ("Lophyrus Abietis", True),   # Synonym of Neodiprion lecontei
        19: ("Lophyrus fulva", True),     # Synonym of Neodiprion edulicolus
        20: ("Xyela ferruginea", False), # Accepted name
        21: ("Xyela aenea", True),        # Synonym of Xyela minor
        22: ("Tremex columba", False)   # Accepted name
    }

    synonym_indices = []
    # Iterate through the species in ascending order of their index
    for index in sorted(species_status.keys()):
        is_synonym = species_status[index][1]
        if is_synonym:
            # Add the index to our list if it's a synonym
            synonym_indices.append(str(index))

    # Join the numbers into a comma-separated string for the final output
    final_output = ",".join(synonym_indices)
    print(final_output)

find_synonym_indices()
<<<2,4,5,6,7,8,9,10,11,12,14,15,17,18,19,21>>>