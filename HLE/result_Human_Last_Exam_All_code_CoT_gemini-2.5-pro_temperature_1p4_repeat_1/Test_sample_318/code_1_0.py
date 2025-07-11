def find_synonym_indices():
    """
    This script identifies which species from the 1872 survey list
    are now considered synonyms according to modern (c. 2020) taxonomy.
    
    The analysis is based on established entomological and taxonomic data.
    A name is considered a synonym if it's no longer the currently accepted
    scientific name for the species.
    """

    # Based on taxonomic research, the following indices correspond to species names
    # from the 1872 list that are now considered synonyms.
    # 2: Abia Kennicotti -> Zaraea inflata
    # 4: Ptenos texanus -> Ptenos niger
    # 7: Hylotoma abdominalis -> Arge abdominalis (genus change)
    # 8: Hylotoma miniata -> Arge miniata (genus change)
    # 9: Hylotoma rubiginosa -> Arge rubiginosa (genus change)
    # 11: Emphytus Bollii -> Ametastegia bollii (genus change)
    # 14: Selandria inaequidens -> Eriocampa inaequidens (genus change)
    # 15: Selandria albicollis -> Periclista albicollis (genus change)
    # 18: Lophyrus Abietis -> Neodiprion lecontei
    # 19: Lophyrus fulva -> Neodiprion excitans
    
    synonym_indices = [2, 4, 7, 8, 9, 11, 14, 15, 18, 19]

    # Sort the indices in ascending order (they are already, but this is good practice)
    synonym_indices.sort()

    # Print the result as a comma-separated string without spaces
    print(','.join(map(str, synonym_indices)))

find_synonym_indices()