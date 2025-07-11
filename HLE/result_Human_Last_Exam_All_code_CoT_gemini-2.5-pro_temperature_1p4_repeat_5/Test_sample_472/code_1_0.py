def scan_latin_line():
    """
    Prints the scansion of a line from Plautus.
    L = Long Syllable
    S = Short Syllable
    Feet are separated by spaces.
    """
    line = "et tibi bene esse soli quom sibi sit male"
    scansion = "L S  S S S  L S  L L  L S  S L  S S"
    
    # Print the original line
    print(f"Line: {line}")
    
    # Print the scansion with feet clearly represented
    print(f"Scansion: {scansion}")
    
    # To fulfill the final output format requirement, we print the numbers in the final equation.
    # The 'equation' here is the scansion itself.
    print("\nFinal Equation:")
    
    syllables_by_foot = scansion.split('  ')
    equation_parts = []
    
    word_syllables = [
        "et", "ti", "bi", "be", "ne", "es", "se", "so", "li", 
        "quom", "si", "bi", "sit", "ma", "le"
    ]
    
    # This is a representation to show which parts of the line correspond to the scansion.
    # It's a symbolic equation.
    
    # Foot 1: et ti -> L S
    equation_parts.append(f"{word_syllables[0]} {word_syllables[1]} -> {syllables_by_foot[0]}")
    # Foot 2: bi be ne -> S S S
    equation_parts.append(f"{word_syllables[2]} {word_syllables[3]} {word_syllables[4]} -> {syllables_by_foot[1]}")
    # Foot 3: es se -> L S
    equation_parts.append(f"{word_syllables[5]} {word_syllables[6]} -> {syllables_by_foot[2]}")
    # Foot 4: so li -> L L
    equation_parts.append(f"{word_syllables[7]} {word_syllables[8]} -> {syllables_by_foot[3]}")
    # Foot 5: quom si -> L S
    equation_parts.append(f"{word_syllables[9]} {word_syllables[10]} -> {syllables_by_foot[4]}")
    # Foot 6: bi sit -> S L
    equation_parts.append(f"{word_syllables[11]} {word_syllables[12]} -> {syllables_by_foot[5]}")
    # Foot 7: ma le -> S S
    equation_parts.append(f"{word_syllables[13]} {word_syllables[14]} -> {syllables_by_foot[6]}")

    final_equation_string = " | ".join(equation_parts)
    
    # Re-printing the final scansion clearly as the 'equation'
    print(scansion.replace("  ", " | "))


scan_latin_line()