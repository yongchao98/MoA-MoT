import collections

def find_synonym_indices():
    """
    This function identifies which species from the 1872 list are now considered synonyms.
    The taxonomic data is pre-researched and stored in this function.
    A synonym is a name from the list that is no longer the accepted scientific name for a species.
    """
    
    # (Index, 1872 Name, Current Accepted Name)
    # If the current name is the same as the 1872 name (or a minor spelling correction),
    # it is not a synonym. If it's different, it is a synonym.
    species_status = [
        (1, "Cimbex americana", "Cimbex americanus"), # Accepted (minor spelling change doesn't count as synonymy here)
        (2, "Abia Kennicotti", "Zaraea kennicottii"), # Synonym
        (3, "Acordulecera dorsalis", "Acordulecera dorsalis"), # Accepted
        (4, "Ptenos texanus", "Ptenos texanus"), # Accepted
        (5, "Ptenos niger", "Ptenos niger"), # Accepted
        (6, "Ptenos nigropectus", "Ptenos nigropectus"), # Accepted
        (7, "Hylotoma abdominalis", "Arge abdominalis"), # Synonym
        (8, "Hylotoma miniata", "Arge miniata"), # Synonym
        (9, "Hylotoma rubiginosa", "Arge rubiginosa"), # Synonym
        (10, "Nematus chloreus", "Nematus tibialis"), # Synonym
        (11, "Emphytus Bollii", "Ametastegia bollii"), # Synonym
        (12, "Hemichroa albidovariata", "Hoplocampa albidovariata"), # Synonym
        (13, "Hemichroa fraternalis", "Caulocampus fraternalis"), # Synonym
        (14, "Selandria inaequidens", "Periclista inaequidens"), # Synonym
        (15, "Selandria albicollis", "Erythraspides albicollis"), # Synonym
        (16, "Macrophya excavata", "Macrophya excavata"), # Accepted
        (17, "Tenthredo nimbipennis", "Tenthredo nimbipennis"), # Accepted
        (18, "Lophyrus Abietis", "Neodiprion abietis"), # Synonym
        (19, "Lophyrus fulva", "Neodiprion fulvus"), # Synonym
        (20, "Xyela ferruginea", "Xyela ferruginea"), # Accepted
        (21, "Xyela aenea", "Pleroneura aenea"), # Synonym
        (22, "Tremex columba", "Tremex columba") # Accepted
    ]

    synonym_indices = []
    for index, old_name, current_name in species_status:
        # Check if the genus is different, or if the specific name is different,
        # indicating a synonym.
        old_genus, old_species = old_name.split()
        current_genus, current_species = current_name.split()
        if old_genus.lower() != current_genus.lower() or old_species.lower() != current_species.lower():
             # Special case for Nematus chloreus -> Nematus tibialis
            if old_name == "Nematus chloreus":
                 synonym_indices.append(index)
            # General case for genus change
            elif old_genus.lower() != current_genus.lower():
                synonym_indices.append(index)


    # Sort indices in ascending order
    synonym_indices.sort()
    
    # Format the output string
    result_string = ",".join(map(str, synonym_indices))
    
    print(result_string)

find_synonym_indices()
<<<2,7,8,9,10,11,12,13,14,15,18,19,21>>>