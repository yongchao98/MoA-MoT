def scan_latin_line():
    """
    Prints the scansion of a line from Plautus.
    L = Long syllable
    S = Short syllable
    """
    line_by_feet = "et ti | bi be | ne es | se so | li quom | si bi | sit ma | le"
    scansion = "L S  | S S  | S L  | S L  | L L   | S S   | L S   | S"
    
    print("Scansion of the line: 'et tibi bene esse soli quom sibi sit male'")
    print("-" * 65)
    print("Syllables grouped by feet:")
    print(line_by_feet)
    print("\nMetrical pattern (L = Long, S = Short):")
    print(scansion)

scan_latin_line()
<<<L S | S S | S L | S L | L L | S S | L S | S>>>