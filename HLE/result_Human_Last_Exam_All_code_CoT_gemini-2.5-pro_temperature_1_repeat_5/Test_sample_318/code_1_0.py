import collections

def find_synonym_indices():
    """
    This function identifies which species from the 1872 list are now considered synonyms.
    The taxonomic data is based on modern databases as of 2020+.
    """
    
    # Step 1: List the species names from the 1872 survey.
    # Note: Minor variations in spelling or capitalization are normalized for consistency.
    species_list_1872 = [
        "Cimbex americana, var. Ulmi", # Accepted: Cimbex americana
        "Abia Kennicotti", # Accepted: Abia kennicottii
        "Acordulecera dorsalis", # Accepted
        "Ptenos texanus", # Synonym
        "Ptenos niger", # Synonym
        "Ptenos nigropectus", # Synonym
        "Hylotoma abdominalis", # Synonym
        "Hylotoma miniata", # Synonym
        "Hylotoma rubiginosa", # Synonym
        "Nematus chloreus", # Synonym
        "Emphytus Bollii", # Synonym
        "Hemichroa albidovariata", # Accepted
        "Hemichroa fraternalis", # Synonym
        "Selandria inaequidens", # Synonym
        "Selandria albicollis", # Synonym
        "Macrophya excavata", # Accepted
        "Tenthredo nimbipennis", # Synonym
        "Lophyrus Abietis", # Synonym
        "Lophyrus fulva", # Synonym
        "Xyela ferruginea", # Accepted
        "Xyela aenea", # Synonym
        "Tremex columba" # Accepted
    ]

    # Step 2 & 3: Define the set of names that are now considered synonyms.
    # This data is based on external taxonomic research.
    # Mapping: {1872 Name -> Modern Accepted Name}
    synonym_map = {
        "Ptenos texanus": "Arge texana",
        "Ptenos niger": "Arge nigra",
        "Ptenos nigropectus": "Arge nigropectus",
        "Hylotoma abdominalis": "Arge abdominalis",
        "Hylotoma miniata": "Arge miniata",
        "Hylotoma rubiginosa": "Arge rubiginosa",
        "Nematus chloreus": "Pontania chlorea",
        "Emphytus Bollii": "Ametastegia bollii",
        "Hemichroa fraternalis": "Hoplocampa fraternalis",
        "Selandria inaequidens": "Caliroa inaequidens",
        "Selandria albicollis": "Erythraspides albicollis",
        "Tenthredo nimbipennis": "Tenthredo anomocerus",
        "Lophyrus Abietis": "Neodiprion sertifer",
        "Lophyrus fulva": "Neodiprion fulvus",
        "Xyela aenea": "Xyela alpigena"
    }

    synonym_names = set(synonym_map.keys())
    
    # Step 4: Find the 1-based indices of the synonyms.
    synonym_indices = []
    for i, name in enumerate(species_list_1872):
        # Extract the core species name for matching
        core_name = name.split(',')[0].strip()
        if core_name in synonym_names:
            synonym_indices.append(i + 1)
            
    # Step 5: Format and print the result.
    result = ",".join(map(str, sorted(synonym_indices)))
    print(result)

find_synonym_indices()