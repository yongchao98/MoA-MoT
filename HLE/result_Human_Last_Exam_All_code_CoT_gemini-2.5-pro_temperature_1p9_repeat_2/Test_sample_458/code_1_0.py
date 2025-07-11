def find_migratory_gliders():
    """
    Identifies dragonfly species from a predefined list that are known
    for a long-distance, soaring/gliding migratory ecology, a trait often
    associated with adaptations like a reduced pterostigma.
    """
    
    # The full list of species provided to the user.
    taxa = {
        1: "Didymops transversa",      # Stream Cruiser - Percher
        2: "Urothemis edwarsi",        # Tropical, strong flier but not a primary glider
        3: "Macrodiplax balteata",     # Marl Pennant - Migratory "glider"
        4: "Pantala flavescens",       # Globe Skimmer - Archetypal migratory "glider"
        5: "Orthetrum cancellatum",    # Black-tailed Skimmer - Percher
        6: "Libelulla quadrimaculata", # Four-spotted Chaser - Migratory "flier/sallier", not "glider"
        7: "Libelulla pulchela",       # Twelve-spotted Skimmer - Percher
        8: "Sympetrum corruptum",      # Variegated Meadowhawk - Migratory, but not in the "glider" guild
        9: "Celithemis elisa",         # Calico Pennant - Percher
        10: "Tholymis tillarga"        # Coral-tailed Cloudwing - Migratory "glider"
    }
    
    # Based on ecological research, species known as long-distance "gliders"
    # are expected to have adaptations for soaring flight, which may include
    # a relatively reduced pterostigma.
    # The indices of these species are 3, 4, and 10.
    reduced_pterostigmata_indices = [3, 4, 10]
    
    # Convert the list of integer indices to a list of strings
    indices_as_strings = [str(index) for index in reduced_pterostigmata_indices]
    
    # Join the string representations with a comma
    result = ",".join(indices_as_strings)
    
    print(result)

find_migratory_gliders()
<<<3,4,10>>>