def find_photosynthetic_taxa():
    """
    Identifies and prints the indices of taxa that perform photochemical synthesis.

    The list of booleans corresponds to the following taxa:
    1) Acanthella cavernosa (False)
    2) Gloeochaete wittrockiana (True)
    3) Homo sapiens (True)
    4) Riftia pachyptila (False)
    5) Halapricum salinum (True)
    6) Aphanothece castagnei (True)
    7) Baileya pleniradiata (True)
    8) Acanthella pulchra (False)
    9) Ectothiorhodosinus mongolicus (True)
    10) Chlorobaculum tepidum (True)
    11) Stygichthys typhlops (False)
    12) Gemmatimonas phototrophica (True)
    13) Myonera garretti (False)
    """
    # Boolean list representing whether each taxon performs photochemical synthesis
    performs_synthesis = [
        False, True, True, False, True, True, True, False, True, True, False, True, False
    ]

    # Find the 1-based indices for the True values
    result_indices = []
    for i, does_synthesize in enumerate(performs_synthesis):
        if does_synthesize:
            # Add 1 to convert 0-based index to 1-based index
            result_indices.append(str(i + 1))

    # Join the indices with a comma and print
    final_answer = ",".join(result_indices)
    print(final_answer)

find_photosynthetic_taxa()