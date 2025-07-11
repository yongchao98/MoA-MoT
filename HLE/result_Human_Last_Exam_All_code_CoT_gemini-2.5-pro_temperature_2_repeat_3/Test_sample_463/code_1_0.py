def solve():
    """
    This function identifies the indices of species from a given list
    that perform photochemical synthesis.
    """
    # List of species provided in the problem
    species_list = [
        "Acanthella cavernosa",        # 1. No (Sponge, animal)
        "Gloeochaete wittrockiana",    # 2. Yes (Alga, photosynthetic)
        "Homo sapiens",                # 3. No (Human, animal)
        "Riftia pachyptila",           # 4. No (Tube worm, relies on chemosynthesis)
        "Halapricum salinum",          # 5. Yes (Archaea, phototrophic)
        "Aphanothece castagnei",       # 6. Yes (Cyanobacteria, photosynthetic)
        "Baileya pleniradiata",        # 7. Yes (Plant, photosynthetic)
        "Acanthella pulchra",          # 8. No (Sponge, animal)
        "Ectothiorhodosinus mongolicus", # 9. Yes (Purple sulfur bacterium, photosynthetic)
        "Chlorobaculum tepidum",       # 10. Yes (Green sulfur bacterium, photosynthetic)
        "Stygichthys typhlops",        # 11. No (Cavefish, animal)
        "Gemmatimonas phototrophica",  # 12. Yes (Bacterium, photosynthetic)
        "Myonera garretti"             # 13. No (Bivalve, animal)
    ]

    # Indices of the species that perform photochemical synthesis
    photochemical_synthesis_indices = [2, 5, 6, 7, 9, 10, 12]

    # Convert the list of numbers to a list of strings
    indices_as_strings = [str(index) for index in photochemical_synthesis_indices]

    # Join the strings with a comma
    result = ",".join(indices_as_strings)

    # Print the final result
    print(result)

solve()