def solve_task():
    """
    This function identifies and prints the indices of species that undergo
    any process of photochemical synthesis as part of their normal metabolic processes,
    ignoring symbiotic relationships.

    The analysis for each species is as follows:
    1) Acanthella cavernosa: No (Heterotrophic animal)
    2) Gloeochaete wittrockiana: Yes (Photosynthetic alga)
    3) Homo sapiens: Yes (Photochemical synthesis of Vitamin D)
    4) Riftia pachyptila: No (Relies on chemosynthesis, no light)
    5) Halapricum salinum: Yes (Phototrophic archaeon, likely uses bacteriorhodopsin)
    6) Aphanothece castagnei: Yes (Photosynthetic cyanobacterium)
    7) Baileya pleniradiata: Yes (Photosynthetic plant)
    8) Acanthella pulchra: No (Heterotrophic animal)
    9) Ectothiorhodosinus mongolicus: Yes (Anoxygenic photosynthetic bacterium)
    10) Chlorobaculum tepidum: Yes (Anoxygenic photosynthetic bacterium)
    11) Stygichthys typhlops: No (Heterotrophic cavefish, no light)
    12) Gemmatimonas phototrophica: Yes (Anoxygenic photosynthetic bacterium)
    13) Myonera garretti: No (Heterotrophic bivalve)
    """

    # List of indices of the species that fit the criteria.
    photochemical_indices = [2, 3, 5, 6, 7, 9, 10, 12]

    # To fulfill the requirement "output each number in the final equation",
    # we will convert the list of numbers into a comma-separated string.
    # First, convert each integer in the list to a string.
    indices_as_strings = [str(index) for index in photochemical_indices]

    # Then, join these strings with a comma.
    final_equation = ",".join(indices_as_strings)

    # Print the final result.
    print(final_equation)

solve_task()