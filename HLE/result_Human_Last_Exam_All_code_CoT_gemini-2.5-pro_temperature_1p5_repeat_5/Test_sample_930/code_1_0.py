def find_mutualists():
    """
    This function identifies the indices of mutualists of Asclepias fascicularis
    from the provided list and prints them.

    - Mutualists are organisms in a relationship where both benefit.
    - Adult pollinators (bees, butterflies, wasps) and defenders (ants) are mutualists.
      They get nectar, and the plant gets pollinated or protected.
    - Herbivorous larvae (e.g., monarch caterpillars) are not mutualists.
    - Larvae that do not interact with the plant (e.g., bee/wasp larvae in a nest) are not mutualists.
    """

    # Indices of organisms identified as mutualists (all the adults)
    mutualist_indices = [1, 2, 3, 4, 5, 6]

    # Convert the list of integers to a list of strings
    mutualist_indices_str = [str(i) for i in mutualist_indices]

    # Join the strings with a comma and print the result
    final_answer = ",".join(mutualist_indices_str)
    print(final_answer)

find_mutualists()