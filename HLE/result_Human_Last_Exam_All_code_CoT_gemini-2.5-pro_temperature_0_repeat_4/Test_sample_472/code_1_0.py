def scan_latin_line():
    """
    This function prints the scansion of the Latin line:
    "et tibi bene esse soli quom sibi sit male"
    L = Long syllable
    S = Short syllable
    Feet are separated by two spaces.
    """
    # The scansion is determined by analyzing the meter (iambic senarius),
    # applying elisions, and determining syllable quantity.
    # Feet: (L S) (S L) (L L) (L S) (L L) (S S)
    scansion = "L S  S L  L L  L S  L L  S S"
    print(scansion)

scan_latin_line()