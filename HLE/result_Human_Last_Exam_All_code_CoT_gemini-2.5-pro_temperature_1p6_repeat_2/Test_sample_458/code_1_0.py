def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species expected to have reduced pterostigmata based on ecology.
    
    The key ecological trait associated with reduced pterostigmata is a flight style
    that involves extensive gliding, typical of long-distance migratory species often
    referred to as "gliders" or "wanderers". This script encodes this knowledge
    to filter the provided list.
    """
    
    # A list of dictionaries holding the ecological data for each species.
    species_data = [
        {'index': 1, 'name': 'Didymops transversa', 'ecology': 'Stream cruiser, a territorial patroller.'},
        {'index': 2, 'name': 'Urothemis edwarsi', 'ecology': 'A migrant, but not typically classified as a primary glider.'},
        {'index': 3, 'name': 'Macrodiplax balteata', 'ecology': 'A strong coastal and trans-oceanic migrant known as a "glider".'},
        {'index': 4, 'name': 'Pantala flavescens', 'ecology': 'The archetypal "wandering glider" and the most widespread dragonfly.'},
        {'index': 5, 'name': 'Orthetrum cancellatum', 'ecology': 'A common territorial percher.'},
        {'index': 6, 'name': 'Libelulla quadrimaculata', 'ecology': 'Known for massive migrations, but primarily a "flapper" not a glider.'},
        {'index': 7, 'name': 'Libelulla pulchela', 'ecology': 'A territorial perching species.'},
        {'index': 8, 'name': 'Sympetrum corruptum', 'ecology': 'A well-known migrant, but not specialized as a glider.'},
        {'index': 9, 'name': 'Celithemis elisa', 'ecology': 'A small percher, not a migrant.'},
        {'index': 10, 'name': 'Tholymis tillarga', 'ecology': 'A widespread crepuscular migrant and "glider", often found with Pantala.'}
    ]
    
    # Keywords that indicate the relevant ecology.
    glider_keywords = ["glider", "wanderer"]
    
    # List to store the indices of the identified species.
    indices_with_reduced_pterostigmata = []
    
    # Iterate through each species and check its ecology.
    for species in species_data:
        # Check if any keyword is present in the ecology description.
        if any(keyword in species['ecology'].lower() for keyword in glider_keywords):
            indices_with_reduced_pterostigmata.append(str(species['index']))
            
    # Format and print the final result.
    if indices_with_reduced_pterostigmata:
        result_string = ",".join(indices_with_reduced_pterostigmata)
        print(result_string)
    else:
        print("none")

find_species_with_reduced_pterostigmata()
<<<3,4,10>>>