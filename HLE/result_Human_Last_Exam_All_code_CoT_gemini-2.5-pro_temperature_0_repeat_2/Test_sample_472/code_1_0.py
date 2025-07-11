def scan_latin_line():
    """
    Prints the scansion of the Latin line.
    L = Long syllable
    S = Short syllable
    Feet are separated by a space.
    """
    # The line is an iambic senarius: et tibi bene esse soli quom sibi sit male
    # Scansion breakdown:
    # Foot 1: et ti bi (L S S) - Dactyl
    # Foot 2: be ne 'sse (S S S) - Tribrach (prodelision of esse)
    # Foot 3: so li (L L) - Spondee
    # Foot 4: quom si (L S) - Trochee
    # Foot 5: bi sit (S L) - Iamb
    # Foot 6: ma le (S S) - Pyrrhic
    
    foot1 = "L S S"
    foot2 = "S S S"
    foot3 = "L L"
    foot4 = "L S"
    foot5 = "S L"
    foot6 = "S S"
    
    # Joining feet with a space as requested
    scansion = " ".join([foot1, foot2, foot3, foot4, foot5, foot6])
    
    print(scansion)

scan_latin_line()