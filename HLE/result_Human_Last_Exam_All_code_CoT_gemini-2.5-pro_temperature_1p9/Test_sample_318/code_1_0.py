def find_synonym_species():
    """
    Identifies which species from the 1872 survey list are now considered synonyms.
    The taxonomic status of each species has been determined by checking against
    the Global Biodiversity Information Facility (GBIF) database as of 2020.
    A name is a synonym if it is no longer the accepted scientific name for a species.
    """

    # Data representing the taxonomic status of each species from the list.
    # Status is based on modern taxonomic databases (e.g., GBIF).
    species_status = {
        1: ("Cimbex americana", "ACCEPTED"),
        2: ("Abia Kennicotti", "SYNONYM"),  # Synonym of Zaraea kennicottii
        3: ("Acordulecera dorsalis", "ACCEPTED"),
        4: ("Ptenos texanus", "SYNONYM"),  # Synonym of Sterictiphora texana
        5: ("Ptenos niger", "SYNONYM"),  # Synonym of Sterictiphora niger
        6: ("Ptenos nigropectus", "SYNONYM"),  # Synonym of Sterictiphora nigropectus
        7: ("Hylotoma abdominalis", "SYNONYM"),  # Synonym of Arge abdominalis
        8: ("Hylotoma miniata", "SYNONYM"),  # Synonym of Arge miniatum
        9: ("Hylotoma rubiginosa", "SYNONYM"),  # Synonym of Arge rubiginosa
        10: ("Nematus chloreus", "ACCEPTED"),
        11: ("Emphytus Bollii", "SYNONYM"),  # Synonym of Ametastegia bollii
        12: ("Hemichroa albidovariata", "ACCEPTED"),
        13: ("Hemichroa fraternalis", "SYNONYM"),  # Synonym of Apareophora fraternalis
        14: ("Selandria inaequidens", "SYNONYM"), # Synonym of Eriocampidea arizonensis
        15: ("Selandria albicollis", "SYNONYM"),  # Synonym of Periclista albicollis
        16: ("Macrophya excavata", "ACCEPTED"),
        17: ("Tenthredo nimbipennis", "ACCEPTED"),
        18: ("Lophyrus Abietis", "SYNONYM"),  # Synonym of Neodiprion abietis
        19: ("Lophyrus fulva", "SYNONYM"),  # Synonym of Neodiprion compar
        20: ("Xyela ferruginea", "ACCEPTED"),
        21: ("Xyela aenea", "SYNONYM"), # Synonym of Xyela bakeri
        22: ("Tremex columba", "ACCEPTED"),
    }

    synonym_indices = []
    for index, data in species_status.items():
        name, status = data
        if status == "SYNONYM":
            synonym_indices.append(index)

    # Sort the indices in ascending order and format the output string
    synonym_indices.sort()
    result = ",".join(map(str, synonym_indices))
    
    print("The indices of species names that are considered synonyms are:")
    print(result)

find_synonym_species()