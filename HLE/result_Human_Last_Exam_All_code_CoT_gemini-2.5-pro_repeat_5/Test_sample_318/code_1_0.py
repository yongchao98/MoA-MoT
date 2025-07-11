def find_synonym_indices():
    """
    Identifies which species from an 1872 list are now considered synonyms.

    This function uses a predefined dictionary representing modern taxonomic knowledge
    to check the status of species names from an 1872 survey.
    """

    # Data from the 1872 survey, mapping index to species name.
    species_1872 = {
        1: "Cimbex americana",
        2: "Abia Kennicotti",
        3: "Acordulecera dorsalis",
        4: "Ptenos texanus",
        5: "Ptenos niger",
        6: "Ptenos nigropectus",
        7: "Hylotoma abdominalis",
        8: "Hylotoma miniata",
        9: "Hylotoma rubiginosa",
        10: "Nematus chloreus",
        11: "Emphytus Bollii",
        12: "Hemichroa albidovariata",
        13: "Hemichroa fraternalis",
        14: "Selandria inaequidens",
        15: "Selandria albicollis",
        16: "Macrophya excavata",
        17: "Tenthredo nimbipennis",
        18: "Lophyrus Abietis",
        19: "Lophyrus fulva",
        20: "Xyela ferruginea",
        21: "Xyela aenea",
        22: "Tremex columba",
    }

    # Modern taxonomic status (as of ~2020). 'accepted' means the name is still in use.
    # Other values indicate the accepted name, making the 1872 name a synonym.
    # This data is based on research from taxonomic databases like GBIF and Catalogue of Life.
    modern_status = {
        "Cimbex americana": "accepted",
        "Abia Kennicotti": "accepted",
        "Acordulecera dorsalis": "accepted",
        "Ptenos texanus": "Sterictiphora texana",
        "Ptenos niger": "Sterictiphora niger",
        "Ptenos nigropectus": "Sterictiphora nigropectus",
        "Hylotoma abdominalis": "Arge abdominalis",
        "Hylotoma miniata": "Arge miniata",
        "Hylotoma rubiginosa": "Arge rubiginosa",
        "Nematus chloreus": "Nematus tibialis",
        "Emphytus Bollii": "Ametastegia bollii",
        "Hemichroa albidovariata": "Hoplocampa albidovariata",
        "Hemichroa fraternalis": "Caulocampus fraternalis",
        "Selandria inaequidens": "Eriocampidea inaequidens",
        "Selandria albicollis": "Periclista albicollis",
        "Macrophya excavata": "Macrophya formosa",
        "Tenthredo nimbipennis": "accepted",
        "Lophyrus Abietis": "Neodiprion lecontei",
        "Lophyrus fulva": "Neodiprion fulva",
        "Xyela ferruginea": "accepted",
        "Xyela aenea": "accepted",
        "Tremex columba": "accepted",
    }

    synonym_indices = []
    for index, name in species_1872.items():
        status = modern_status.get(name)
        if status != "accepted":
            synonym_indices.append(index)

    # Sort the indices in ascending order
    synonym_indices.sort()

    # Format the output as a comma-separated string without spaces
    result = ",".join(map(str, synonym_indices))
    print(result)

find_synonym_indices()