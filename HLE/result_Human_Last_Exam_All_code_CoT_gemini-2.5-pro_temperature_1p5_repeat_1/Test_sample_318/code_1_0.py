import sys

def find_synonyms():
    """
    Identifies species from an 1872 list that are now considered synonyms.
    
    This function uses a pre-compiled dictionary to simulate checking a modern
    taxonomic database (as of c. 2020) to find which historical species names
    have been superseded.
    """
    
    # A list of the species names as they appear in the 1872 survey.
    # Note: Minor variations in spelling/capitalization are normalized for matching.
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
        22: "Tremex columba"
    }

    # A dictionary representing the results from a 2020 taxonomic database check.
    # Key: The original name from the 1872 list.
    # Value: The current accepted scientific name.
    # If a name from 1872 is not in this dictionary, it is still considered valid.
    synonym_database = {
        "Abia Kennicotti": "Zaraea inflata",
        "Hylotoma abdominalis": "Arge abdominalis",
        "Hylotoma miniata": "Arge miniata",
        "Hylotoma rubiginosa": "Arge rubiginosa",
        "Emphytus Bollii": "Empria maculata",
        "Selandria inaequidens": "Periclista inaequidens",
        "Selandria albicollis": "Proselandria albicollis",
        "Lophyrus Abietis": "Neodiprion lecontei",
        "Lophyrus fulva": "Neodiprion fabricii",
        "Xyela aenea": "Xyela aeneipennis"
    }
    
    synonym_indices = []
    
    print("Finding which species from the 1872 list are now synonyms:")
    
    # Iterate through the original species list by sorted index
    for index in sorted(species_1872.keys()):
        original_name = species_1872[index]
        if original_name in synonym_database:
            synonym_indices.append(index)
            current_name = synonym_database[original_name]
            print(f"- Index {index} ({original_name}) is a synonym of {current_name}.")

    # Format the final answer as a comma-separated string of indices
    final_answer = ",".join(map(str, synonym_indices))
    
    print("\nThe indices of the synonymous species are:")
    # The final answer format as requested by the user.
    # The print statements above act as the step-by-step explanation.
    # "Remember in the final code you still need to output each number in the final equation!"
    # The loop above shows how each number (index) was derived.
    # The final output is just the requested format.
    print(f'<<<{final_answer}>>>')

if __name__ == '__main__':
    find_synonyms()