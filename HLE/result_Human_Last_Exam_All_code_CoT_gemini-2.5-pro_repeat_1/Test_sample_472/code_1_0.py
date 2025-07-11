def scan_latin_line():
    """
    Scans the provided Latin line and prints the result.
    The line is an Iambic Senarius from Plautus.
    Line: "et tibi bene esse soli quom sibi sit male"

    Scansion Breakdown:
    - The line is scanned as 13 syllables after two elisions: `tibi bene` -> `tib'` and `bene esse` -> `benesse`.
    - The resulting line for scanning is: "et tib' benesse soli quom sibi sit male"

    Foot 1: et tib' -> L S (Trochee)
    Foot 2: be-nes-se -> S L L (Bacchius, a resolved spondee)
    Foot 3: so-li -> L L (Spondee)
    Foot 4: quom si -> L S (Trochee)
    Foot 5: -bi sit -> S L (Iamb)
    Foot 6: ma-le -> S L (Iamb, due to 'brevis in longo')
    """
    foot1 = "L S"
    foot2 = "S L L"
    foot3 = "L L"
    foot4 = "L S"
    foot5 = "S L"
    foot6 = "S L"
    
    # Print the scansion with a space between each foot
    print(foot1, foot2, foot3, foot4, foot5, foot6)

scan_latin_line()