def scan_latin_line():
    """
    This function provides the scansion for the given Latin line.
    L = Long syllable
    S = Short syllable
    Feet are separated by spaces.
    """
    # The line: "et tibi bene esse soli quom sibi sit male"
    # The scansion is a highly substituted Iambic Senarius.
    # Foot 1: L S (Trochee)
    # Foot 2: S S (Pyrrhic)
    # Foot 3: L S (Iamb)
    # Foot 4: L L (Spondee)
    # Foot 5: L S S (Dactyl)
    # Foot 6: L S S (Dactyl)
    scansion = "L S S S L S L L L S S L S S"
    
    # Adding spaces to separate the feet
    feet = ["L S", "S S", "L S", "L L", "L S S", "L S S"]
    final_scansion = " ".join(feet)
    
    print(final_scansion)

scan_latin_line()