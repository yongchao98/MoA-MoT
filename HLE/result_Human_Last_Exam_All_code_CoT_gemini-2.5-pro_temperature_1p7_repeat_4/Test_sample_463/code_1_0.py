def solve():
    """
    This function identifies and prints the indices of taxa that perform photochemical synthesis.
    The analysis of each taxon is as follows:
    1) Acanthella cavernosa: No (Sponge, heterotroph)
    2) Gloeochaete wittrockiana: Yes (Alga, photosynthesis)
    3) Homo sapiens: Yes (Vitamin D synthesis via UV light)
    4) Riftia pachyptila: No (Tube worm, chemosynthesis)
    5) Halapricum salinum: Yes (Archaea, bacteriorhodopsin-based phototrophy)
    6) Aphanothece castagnei: Yes (Cyanobacteria, photosynthesis)
    7) Baileya pleniradiata: Yes (Plant, photosynthesis)
    8) Acanthella pulchra: No (Sponge, heterotroph)
    9) Ectothiorhodosinus mongolicus: Yes (Purple sulfur bacteria, anoxygenic photosynthesis)
    10) Chlorobaculum tepidum: Yes (Green sulfur bacteria, anoxygenic photosynthesis)
    11) Stygichthys typhlops: No (Cavefish, heterotroph)
    12) Gemmatimonas phototrophica: Yes (Bacteria, anoxygenic photosynthesis)
    13) Myonera garretti: No (Bivalve, heterotroph)
    """
    photochemical_indices = [2, 3, 5, 6, 7, 9, 10, 12]
    # Convert list of integers to a comma-separated string for printing
    result_string = ",".join(map(str, photochemical_indices))
    print(result_string)

solve()