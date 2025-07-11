def find_dragonflies_with_reduced_pterostigmata():
    """
    Identifies dragonflies from a list that are expected to have reduced 
    pterostigmata based on their ecology.

    The primary ecological trait associated with reduced pterostigmata is being a
    highly specialized, long-distance migratory glider. This script encodes this
    biological knowledge to identify the relevant species.
    """
    
    # Data representing known dragonfly ecologies
    species_data = [
        {"index": 1, "name": "Didymops transversa", "trait": "Normal pterostigmata"},
        {"index": 2, "name": "Urothemis edwarsi", "trait": "Normal pterostigmata"},
        {"index": 3, "name": "Macrodiplax balteata", "trait": "Migratory glider, reduced pterostigmata"},
        {"index": 4, "name": "Pantala flavescens", "trait": "Extreme migratory glider, reduced pterostigmata"},
        {"index": 5, "name": "Orthetrum cancellatum", "trait": "Normal pterostigmata"},
        {"index": 6, "name": "Libelulla quadrimaculata", "trait": "Migratory, but normal pterostigmata"},
        {"index": 7, "name": "Libelulla pulchela", "trait": "Prominent pterostigmata"},
        {"index": 8, "name": "Sympetrum corruptum", "trait": "Migratory, but normal pterostigmata"},
        {"index": 9, "name": "Celithemis elisa", "trait": "Prominent pterostigmata"},
        {"index": 10, "name": "Tholymis tillarga", "trait": "Migratory glider, reduced pterostigmata"},
    ]

    print("Identifying species with reduced pterostigmata based on ecology (long-distance gliding):")
    
    target_indices = []
    for species in species_data:
        # Check for the specific ecological trait
        if "reduced pterostigmata" in species["trait"]:
            target_indices.append(str(species["index"]))
            print(f"- Species {species['index']} ({species['name']}): Identified as a migratory glider with reduced pterostigmata.")

    if not target_indices:
        print("none")
    else:
        # The prompt asks to "output each number in the final equation", which we interpret as printing the final list.
        final_answer = ",".join(target_indices)
        print("\nFinal list of indices:")
        print(final_answer)

if __name__ == "__main__":
    find_dragonflies_with_reduced_pterostigmata()