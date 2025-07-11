def scan_latin_line():
    """
    Scans the provided Latin line, an iambic senarius, and prints the result.
    The scansion is based on expert analysis, including elisions and metrical substitutions.
    """
    
    # The line is "et tibi bene esse soli quom sibi sit male".
    # After elision (tib(i) bene and ben(e) esse), the pronounced syllables and their quantities are:
    
    # Syllables with their determined quantities (L=long, S=short)
    syllables = [
        # Word: et
        ('et', 'L'),  # Long by position before 't' of tibi
        # Word: tibi (elided to tib)
        ('tib', 'L'), # Long by position before 'b' of bene
        # Word: bene (elided to ben)
        ('ben', 'L'), # Long by position before 's' of esse
        # Word: esse
        ('es', 'L'),  # Long by position (ss)
        ('se', 'S'),  # Naturally short
        # Word: soli
        ('so', 'L'),  # Naturally long 'o'
        ('li', 'S'),  # Short 'i'
        # Word: quom
        ('quom', 'L'), # Long by position before 's' of sibi
        # Word: sibi
        ('si', 'S'),  # Naturally short
        ('bi', 'S'),  # Naturally short (via iambic shortening)
        # Word: sit
        ('sit', 'L'), # Long by position before 'm' of male
        # Word: male
        ('ma', 'S'),  # Naturally short
        ('le', 'L')   # Final syllable of the line, treated as long
    ]
    
    # The meter is Iambic Senarius (6 feet).
    # The 13 syllables are grouped into feet with metrical substitutions.
    # Feet: Spondee | Spondee | Iamb | Iamb | Anapest | Iamb
    foot_structure = [2, 2, 2, 2, 3, 2]
    
    # Print the scansion with original words and foot breaks
    print("Scansion of the line:")
    print("\"et tibi bene esse soli quom sibi sit male\"")
    print("-" * 30)

    syllable_pointer = 0
    full_scansion_string = ""
    # We still need to print the full equation as requested
    final_equation_parts = []
    
    # 1. First foot (Spondee)
    foot1 = syllables[0][1] + " " + syllables[1][1]
    final_equation_parts.append(foot1)
    
    # 2. Second foot (Spondee)
    foot2 = syllables[2][1] + " " + syllables[3][1]
    final_equation_parts.append(foot2)
    
    # 3. Third foot (Iamb)
    foot3 = syllables[4][1] + " " + syllables[5][1]
    final_equation_parts.append(foot3)
    
    # 4. Fourth foot (Iamb)
    foot4 = syllables[6][1] + " " + syllables[7][1]
    final_equation_parts.append(foot4)

    # 5. Fifth foot (Anapest)
    foot5 = syllables[8][1] + " " + syllables[9][1] + " " + syllables[10][1]
    final_equation_parts.append(foot5)

    # 6. Sixth foot (Iamb)
    foot6 = syllables[11][1] + " " + syllables[12][1]
    final_equation_parts.append(foot6)

    print(" | ".join(final_equation_parts))

scan_latin_line()