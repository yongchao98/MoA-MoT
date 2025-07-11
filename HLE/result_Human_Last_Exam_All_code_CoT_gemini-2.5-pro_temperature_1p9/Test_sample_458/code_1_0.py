import sys

def solve_dragonfly_puzzle():
    """
    Identifies dragonfly species with reduced pterostigmata based on their ecology.

    The primary ecological driver for reduced pterostigmata in dragonflies is
    a lifestyle that involves long-distance migration and extensive gliding.
    Species adapted for this 'glider' lifestyle have evolved specific aerodynamic
    traits, including broad hindwings and smaller pterostigmata, for maximum
    efficiency during sustained flight.

    We will evaluate the provided list of species against this ecological profile.
    """
    
    species_data = {
        1: ("Didymops transversa", "Stream Cruiser - A stream patroller, not a long-distance glider."),
        2: ("Urothemis edwarsi", "A tropical species, not primarily known for the kind of specialized gliding flight that leads to reduced pterostigmata."),
        3: ("Macrodiplax balteata", "Marl Pennant - A known long-distance coastal migrant with a gliding flight style, similar to Pantala."),
        4: ("Pantala flavescens", "Wandering Glider - The archetypal example of a global migrant with highly specialized wings for gliding, including reduced pterostigmata."),
        5: ("Orthetrum cancellatum", "Black-tailed Skimmer - Primarily a 'percher' species, does not undertake extensive gliding migrations."),
        6: ("Libelulla quadrimaculata", "Four-spotted Chaser - A known migrant, but not a specialized glider like Pantala; retains more typical libellulid wing features."),
        7: ("Libelulla pulchela", "Twelve-spotted Skimmer - A common percher and territorial patroller, not a long-distance glider."),
        8: ("Sympetrum corruptum", "Variegated Meadowhawk - A migratory species, but not specialized for gliding to the same extent as Pantala or Tholymis."),
        9: ("Celithemis elisa", "Calico Pennant - A 'percher' typically found near ponds, not a long-distance migrant."),
        10: ("Tholymis tillarga", "Coral-tailed Cloudwing - A well-known crepuscular migrant that uses a gliding flight style over long distances.")
    }

    print("Analysis of dragonfly species based on flight ecology:")
    
    glider_indices = []
    for index, (name, reason) in species_data.items():
        is_glider = "glider" in reason or "gliding" in reason
        if is_glider:
            glider_indices.append(index)
            print(f"- Species {index} ({name}): Identified as a long-distance migrant and glider. Expected to have reduced pterostigmata.")
        else:
            print(f"- Species {index} ({name}): Not a specialized glider. Not expected to have reduced pterostigmata.")

    # Sort indices for consistent output
    glider_indices.sort()

    # Format the final answer as a comma-separated string
    final_answer = ",".join(map(str, glider_indices))
    
    print("\nConclusion: The species expected to have reduced pterostigmata are the specialized long-distance gliders.")
    print("The indices of these species are:")
    print(final_answer)

solve_dragonfly_puzzle()