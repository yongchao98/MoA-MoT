def solve_photochemical_synthesis():
    """
    This function identifies species from a predefined list that perform photochemical synthesis.

    The logic is based on biological classification:
    - 2: Gloeochaete wittrockiana is a green alga, which performs oxygenic photosynthesis.
    - 5: Halapricum salinum is a haloarchaeon that uses bacteriorhodopsin for phototrophy, a non-chlorophyll-based photochemical process.
    - 6: Aphanothece castagnei is a cyanobacterium, which performs oxygenic photosynthesis.
    - 7: Baileya pleniradiata is a desert marigold (plant), which performs oxygenic photosynthesis.
    - 9: Ectothiorhodosinus mongolicus is a purple sulfur bacterium, which performs anoxygenic photosynthesis.
    - 10: Chlorobaculum tepidum is a green sulfur bacterium, which performs anoxygenic photosynthesis.
    - 12: Gemmatimonas phototrophica is a bacterium known to perform anoxygenic photosynthesis.

    Other species are excluded:
    - 1, 8 (Sponges), 3 (Human), 4 (Tube worm), 11 (Fish), 13 (Bivalve) are animals and do not perform this process intrinsically.
    """

    # Indices of the species that perform photochemical synthesis
    photosynthetic_indices = [2, 5, 6, 7, 9, 10, 12]

    # Convert the list of integers to a list of strings
    indices_as_strings = [str(i) for i in photosynthetic_indices]

    # Join the strings with a comma to create the final output format
    result = ",".join(indices_as_strings)

    # Print the final result
    print(result)

solve_photochemical_synthesis()