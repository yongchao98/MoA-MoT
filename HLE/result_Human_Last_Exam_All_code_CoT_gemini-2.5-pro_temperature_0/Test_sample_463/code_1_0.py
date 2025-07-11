def solve_task():
    """
    This function identifies and prints the indices of taxa that undergo photochemical synthesis.

    The analysis for each taxon is as follows:
    1) Acanthella cavernosa: No. Animal, heterotroph.
    2) Gloeochaete wittrockiana: Yes. Alga, performs oxygenic photosynthesis.
    3) Homo sapiens: Yes. Animal, performs photochemical synthesis of Vitamin D.
    4) Riftia pachyptila: No. Animal, relies on chemosynthesis, not photosynthesis.
    5) Halapricum salinum: Yes. Archaea, performs bacteriorhodopsin-based phototrophy (ATP synthesis).
    6) Aphanothece castagnei: Yes. Cyanobacteria, performs oxygenic photosynthesis.
    7) Baileya pleniradiata: Yes. Plant, performs oxygenic photosynthesis.
    8) Acanthella pulchra: No. Animal, heterotroph.
    9) Ectothiorhodosinus mongolicus: Yes. Purple sulfur bacteria, performs anoxygenic photosynthesis.
    10) Chlorobaculum tepidum: Yes. Green sulfur bacteria, performs anoxygenic photosynthesis.
    11) Stygichthys typhlops: No. Animal (cavefish), lives in darkness.
    12) Gemmatimonas phototrophica: Yes. Bacterium, performs anoxygenic photosynthesis.
    13) Myonera garretti: No. Animal (bivalve), heterotroph.
    """
    
    # List of indices for the species that perform photochemical synthesis.
    qualifying_indices = [2, 3, 5, 6, 7, 9, 10, 12]
    
    # Convert the list of numbers to a list of strings
    string_indices = [str(i) for i in qualifying_indices]
    
    # Join the strings with a comma to form the final output string
    result = ",".join(string_indices)
    
    # Print the final result
    print(result)

solve_task()