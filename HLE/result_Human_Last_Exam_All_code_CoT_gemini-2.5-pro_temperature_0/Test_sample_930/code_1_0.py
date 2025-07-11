def find_mutualists():
    """
    This function identifies and prints the indices of organisms that are 
    mutualists of Asclepias fascicularis from a predefined list.

    A mutualistic relationship is one where both organisms benefit. In this case,
    it refers to pollinators or protectors of the plant.

    - Adults (1-6): All listed adults (butterflies, bees, wasps, ants) visit the
      plant for nectar and act as pollinators or protectors. This is a mutualistic relationship.
    - Larvae (7-12): The monarch larva (7) is an herbivore. The other larvae (8-12)
      do not interact with the plant. None of these are mutualistic relationships.
    """
    
    # Indices of the organisms that are mutualists
    mutualist_indices = [1, 2, 3, 4, 5, 6]
    
    # Convert the list of numbers to a comma-separated string for printing
    output = ",".join(map(str, mutualist_indices))
    
    print(output)

find_mutualists()