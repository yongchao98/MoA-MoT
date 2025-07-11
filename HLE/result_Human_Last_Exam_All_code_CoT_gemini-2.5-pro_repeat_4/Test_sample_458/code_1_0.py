def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species expected to have reduced pterostigmata based on their ecology.

    The key ecological trait associated with a reduced pterostigma is long-distance,
    gliding migration. This adaptation reduces drag and weight at the wingtip,
    improving aerodynamic efficiency for sustained flight.

    This function analyzes a predefined list of species and selects those known
    for this migratory behavior.
    """

    species_data = {
        1: {'name': 'Didymops transversa', 'is_migratory_glider': False},
        2: {'name': 'Urothemis edwarsi', 'is_migratory_glider': False},
        3: {'name': 'Macrodiplax balteata', 'is_migratory_glider': True},
        4: {'name': 'Pantala flavescens', 'is_migratory_glider': True},
        5: {'name': 'Orthetrum cancellatum', 'is_migratory_glider': False},
        6: {'name': 'Libelulla quadrimaculata', 'is_migratory_glider': False},
        7: {'name': 'Libelulla pulchela', 'is_migratory_glider': False},
        8: {'name': 'Sympetrum corruptum', 'is_migratory_glider': False}, # Migratory, but not a specialized glider in the same vein as Pantala.
        9: {'name': 'Celithemis elisa', 'is_migratory_glider': False},
        10: {'name': 'Tholymis tillarga', 'is_migratory_glider': True}
    }

    print("Analysis: Species with reduced pterostigmata are typically long-distance gliders.")
    print("Identifying species classified as 'migratory_glider'...")
    
    reduced_pterostigmata_indices = []
    for index, data in species_data.items():
        if data['is_migratory_glider']:
            reduced_pterostigmata_indices.append(str(index))
            print(f"- Species {index} ({data['name']}): Identified as a long-distance migratory glider.")

    final_answer = ",".join(reduced_pterostigmata_indices)
    
    print("\nFinal list of indices for species expected to have reduced pterostigmata:")
    print(final_answer)

find_species_with_reduced_pterostigmata()