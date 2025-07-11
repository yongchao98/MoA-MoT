def scan_latin_line():
    """
    Prints the scansion of a given Latin line from Catullus.
    L = Long syllable
    S = Short syllable
    """
    line = "et tibi bene esse soli quom sibi sit male"
    
    # The words are grouped by metrical foot based on the analysis.
    feet_words = ["et tibi", "bene esse", "soli", "quom sibi", "sit male"]
    
    # The scansion for each foot.
    feet_scansion = ["L S S", "L S S", "L L", "L S S", "L S S"]

    print(f"Line: {line}\n")
    print("Scansion:")

    # Print the feet with the corresponding words and scansion aligned.
    header = " | ".join(f"{word:<10}" for word in feet_words)
    scansion_line = " | ".join(f"{scansion:<10}" for scansion in feet_scansion)
    
    print(header)
    print(scansion_line)

scan_latin_line()
<<<L S S | L S S | L L | L S S | L S S>>>