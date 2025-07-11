def find_synonym_indices():
    """
    This function identifies and prints the indices of species from the 1872 list
    that are now considered synonyms based on modern taxonomic data.
    
    The research and determination of synonymy have been done prior to writing this code.
    The code serves to format and present the final answer.
    """
    
    # Based on research from taxonomic databases (e.g., GBIF, ITIS),
    # the following indices correspond to species names from the 1872 list
    # that are now considered synonyms of other species.
    # 1. Cimbex americana, var. Ulmi -> Cimbex americana
    # 2. Abia Kennicotti -> Zaraea kennicottii
    # 4. Ptenos texanus -> Sterictiphora texana
    # 5. Ptenos niger -> Sterictiphora niger
    # 6. Ptenos nigropectus -> Sterictiphora nigropectus
    # 7. Hylotoma abdominalis -> Arge abdominalis
    # 8. Hylotoma miniata -> Arge miniata
    # 9. Hylotoma rubiginosa -> Arge rubiginosa
    # 11. Emphytus Bollii -> Ametastegia bollii
    # 13. Hemichroa fraternalis -> Hemichroa albidovariata
    # 14. Selandria inaequidens -> Caliroa inaequidens
    # 15. Selandria albicollis -> Caliroa albicollis
    # 16. Macrophya excavata -> Macrophya formosa
    # 18. Lophyrus Abietis -> Neodiprion lecontei
    # 19. Lophyrus fulva -> Neodiprion fulviceps
    # 21. Xyela aenea (Norton, 1872) -> Xyela aenea (Say, 1824)
    
    synonym_indices = [1, 2, 4, 5, 6, 7, 8, 9, 11, 13, 14, 15, 16, 18, 19, 21]
    
    # Sort the indices in ascending order (they are already sorted)
    synonym_indices.sort()
    
    # Format the list into a comma-separated string without spaces
    output_string = ",".join(map(str, synonym_indices))
    
    # Print the final result
    print(output_string)

find_synonym_indices()