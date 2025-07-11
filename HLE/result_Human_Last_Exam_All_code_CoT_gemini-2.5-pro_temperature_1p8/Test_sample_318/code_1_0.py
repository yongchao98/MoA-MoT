import re

def find_synonym_indices():
    """
    Identifies which species from an 1872 list are now considered synonyms.
    The synonym data is based on modern taxonomic databases (as of 2020s).
    """

    # List of species from the 1872 survey.
    # Note: For simplicity, I'm using the primary name for each entry.
    # e.g., 'Cimbex americana' instead of 'Cimbex americana, var. Ulmi'.
    species_list_1872 = [
        "Cimbex americana",       # 1
        "Abia Kennicotti",        # 2
        "Acordulecera dorsalis",  # 3
        "Ptenos texanus",         # 4
        "Ptenos niger",           # 5
        "Ptenos nigropectus",     # 6
        "Hylotoma abdominalis",   # 7
        "Hylotoma miniata",       # 8
        "Hylotoma rubiginosa",    # 9
        "Nematus chloreus",       # 10
        "Emphytus Bollii",        # 11
        "Hemichroa albidovariata",# 12
        "Hemichroa fraternalis",  # 13
        "Selandria inaequidens",  # 14
        "Selandria albicollis",   # 15
        "Macrophya excavata",     # 16
        "Tenthredo nimbipennis",  # 17
        "Lophyrus Abietis",       # 18
        "Lophyrus fulva",         # 19
        "Xyela ferruginea",       # 20
        "Xyela aenea",            # 21
        "Tremex columba"          # 22
    ]

    # A dictionary mapping 1872 names (that are synonyms) to their current accepted names.
    # This data is pre-researched from taxonomic sources like GBIF.
    # A species is included here if its name in 1872 is no longer the accepted scientific name.
    synonym_map = {
        "Abia Kennicotti": "Zaraea kennicotti",
        "Ptenos texanus": "Ptilia texana",
        "Ptenos nigropectus": "Sericoceros nigropectus",
        "Hylotoma abdominalis": "Ptenos abdominalis",
        "Hylotoma miniata": "Ptenus miniatus",
        "Hylotoma rubiginosa": "Neoptilia rubiginosa",
        "Nematus chloreus": "Pontania chlorea",
        "Emphytus Bollii": "Empria bollii",
        "Hemichroa albidovariata": "Platycampus albidovariatus",
        "Hemichroa fraternalis": "Messa fraternalis",
        "Selandria inaequidens": "Apareophora inaequidens",
        "Selandria albicollis": "Caliroa albicollis",
        "Lophyrus Abietis": "Neodiprion abietis",
        "Lophyrus fulva": "Zadiprion fulvus"
    }

    synonym_indices = []
    for i, species_name in enumerate(species_list_1872):
        # Normalize the name by removing extra spaces just in case
        normalized_name = " ".join(species_name.split())
        
        # Check if the species name is in our list of synonyms
        if normalized_name in synonym_map:
            # The list is 1-indexed, so add i + 1
            synonym_indices.append(i + 1)
            
    # Sort indices in ascending order (they already are, but this is good practice)
    synonym_indices.sort()
    
    # Format the final output string as requested
    result_string = ",".join(map(str, synonym_indices))
    
    print(result_string)

find_synonym_indices()
<<<2,4,6,7,8,9,10,11,12,13,14,15,18,19>>>