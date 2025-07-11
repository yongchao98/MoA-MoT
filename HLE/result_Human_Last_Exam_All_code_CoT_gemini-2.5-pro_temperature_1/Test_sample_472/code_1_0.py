def scan_latin_line():
    """
    This function prints the scansion of the given Latin line.
    The line is an iambic senarius from Plautus.
    L = Long syllable
    S = Short syllable
    Feet are separated by a space.
    """
    
    # The scansion is determined by analyzing elision, syllable weight, and metrical feet.
    # The line is: "et tibi bene esse soli quom sibi sit male"
    # Elision: bene esse -> ben'esse
    # Feet breakdown:
    # 1. et ti | (L S) - Trochee
    # 2. bi es | (S L) - Iamb
    # 3. se so | (L S) - Trochee
    # 4. li quom | (S L) - Iamb
    # 5. si bi sit | (S S L) - Anapest
    # 6. ma le | (S S) - Pyrrhic (brevis in longo)
    
    foot1 = "L S"
    foot2 = "S L"
    foot3 = "L S"
    foot4 = "S L"
    foot5 = "S S L"
    foot6 = "S S"
    
    final_scansion = f"{foot1} {foot2} {foot3} {foot4} {foot5} {foot6}"
    
    print(final_scansion)

scan_latin_line()
<<<L S S L L S S L S S L S S>>>