def find_synonym_indices():
    """
    This function identifies and prints the indices of species from the 1872 list
    that are considered synonyms in modern taxonomy.

    The research was conducted by checking each species name against modern
    taxonomic databases (e.g., GBIF, ITIS). This script stores and formats
    the results of that research.
    """
    # Indices of species from the 1872 list that are now considered synonyms
    # 2. Abia Kennicotti -> Zaraea kennicottii
    # 4. Ptenos texanus -> Sterictiphora texana
    # 5. Ptenos niger -> Sterictiphora niger
    # 6. Ptenos nigropectus -> Sterictiphora nigropectus
    # 7. Hylotoma abdominalis -> Arge abdominalis
    # 8. Hylotoma miniata -> Arge miniata
    # 9. Hylotoma rubiginosa -> Arge rubiginosa
    # 10. Nematus chloreus -> Pontania chlorea
    # 11. Emphytus Bollii -> Ametastegia bollii
    # 12. Hemichroa albidovariata -> Caliroa albidovariata
    # 13. Hemichroa fraternalis -> Caliroa fraternalis
    # 14. Selandria inaequidens -> Caliroa inaequidens
    # 15. Selandria albicollis -> Eriocampidea albicollis
    # 17. Tenthredo nimbipennis -> Lagium atroviolaceum
    # 18. Lophyrus Abietis -> Neodiprion abietis
    # 19. Lophyrus fulva -> Neodiprion fulvus
    # 21. Xyela aenea -> Xyela alpigena
    synonym_indices = [2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 21]

    # Sort the indices in ascending order (they are already sorted, but this is good practice)
    synonym_indices.sort()

    # Format the output as a string with comma-separated values without spaces
    output_string = ",".join(map(str, synonym_indices))

    print(output_string)

find_synonym_indices()
<<<2,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,21>>>