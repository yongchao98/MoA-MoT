def find_mutualists():
    """
    This function identifies the mutualists of Asclepias fascicularis from a predefined list.
    
    A mutualistic relationship is one where both organisms benefit. For a plant like milkweed,
    this primarily involves pollination or defense.

    - Adult butterflies, bees, and wasps are pollinators, feeding on nectar and transferring pollen.
    - Adult ants feed on nectar and can defend the plant from herbivores.
    - Larvae of these species are either herbivores (harming the plant) or do not interact with the plant directly.
    """
    
    # Indices of organisms that are mutualists (pollinators/defenders)
    mutualist_indices = [1, 2, 3, 4, 5, 6]
    
    # Convert the list of numbers to a list of strings to use the join method
    result_list = [str(index) for index in mutualist_indices]
    
    # Join the list elements into a single string, separated by commas
    final_answer = ",".join(result_list)
    
    print(final_answer)

find_mutualists()