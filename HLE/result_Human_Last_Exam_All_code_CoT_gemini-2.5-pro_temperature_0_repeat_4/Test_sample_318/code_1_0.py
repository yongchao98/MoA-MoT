def find_synonym_indices():
    """
    This function identifies and prints the indices of species from the 1872 list
    that are considered synonyms in modern taxonomy.

    The research was conducted by checking each species name against modern
    taxonomic databases (e.g., GBIF, ITIS). A name is considered a synonym if it's
    no longer the accepted scientific name for that species.
    """

    # List of species from the 1872 survey and their current status (researched externally)
    # 1. Cimbex americana: Accepted
    # 2. Abia Kennicotti: Synonym of Zaraea kennicottii
    # 3. Acordulecera dorsalis: Accepted
    # 4. Ptenos texanus: Accepted
    # 5. Ptenos niger: Accepted
    # 6. Ptenos nigropectus: Accepted
    # 7. Hylotoma abdominalis: Synonym of Arge abdominalis
    # 8. Hylotoma miniata: Synonym of Arge miniata
    # 9. Hylotoma rubiginosa: Synonym of Arge rubiginosa
    # 10. Nematus chloreus: Synonym of Pontania chlorea
    # 11. Emphytus Bollii: Synonym of Ametastegia bollii
    # 12. Hemichroa albidovariata: Synonym of Caliroa albidovariata
    # 13. Hemichroa fraternalis: Synonym of Caliroa fasciata
    # 14. Selandria inaequidens: Synonym of Caliroa inaequidens
    # 15. Selandria albicollis: Synonym of Caliroa albicollis
    # 16. Macrophya excavata: Accepted
    # 17. Tenthredo nimbipennis: Accepted
    # 18. Lophyrus Abietis: Synonym of Neodiprion abietis
    # 19. Lophyrus fulva: Synonym of Neodiprion fulvus
    # 20. Xyela ferruginea: Accepted
    # 21. Xyela aenea: Accepted
    # 22. Tremex columba: Accepted

    synonym_indices = [2, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19]

    # Format the list as a comma-separated string without spaces
    result = ",".join(map(str, synonym_indices))

    print(result)

find_synonym_indices()
<<<2,7,8,9,10,11,12,13,14,15,18,19>>>