def find_photochemical_synthesizers():
    """
    Identifies and prints the indices of taxa that perform photochemical synthesis.

    The function is based on the following analysis:
    - 2) Gloeochaete wittrockiana: An alga that performs oxygenic photosynthesis.
    - 5) Halapricum salinum: An archaeon that performs retinal-based phototrophy (using bacteriorhodopsin to synthesize ATP).
    - 6) Aphanothece castagnei: A cyanobacterium that performs oxygenic photosynthesis.
    - 7) Baileya pleniradiata: A plant that performs oxygenic photosynthesis.
    - 9) Ectothiorhodosinus mongolicus: A purple sulfur bacterium that performs anoxygenic photosynthesis.
    - 10) Chlorobaculum tepidum: A green sulfur bacterium that performs anoxygenic photosynthesis.
    - 12) Gemmatimonas phototrophica: A bacterium that performs anoxygenic photosynthesis.

    The other taxa are non-photosynthetic animals (1, 3, 4, 8, 11, 13).
    """
    phototroph_indices = [2, 5, 6, 7, 9, 10, 12]

    # Convert indices to strings to use the join method
    phototroph_indices_str = [str(i) for i in phototroph_indices]

    # Join the strings with a comma and print the result
    result = ",".join(phototroph_indices_str)
    print(result)

find_photochemical_synthesizers()