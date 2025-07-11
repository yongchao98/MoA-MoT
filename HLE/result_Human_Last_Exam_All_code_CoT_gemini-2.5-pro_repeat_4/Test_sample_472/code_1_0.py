def scan_latin_line():
    """
    This function stores and prints the metrical scansion of the given Latin line.
    L represents a long syllable.
    S represents a short syllable.
    Each group of letters represents a metrical foot.
    """
    # The scansion is determined by applying the rules of Latin prosody.
    # The line is: "et tibi bene esse soli quom sibi sit male"
    # It scans into 6 feet as follows:
    foot1 = "L S S"  # et ti-bi
    foot2 = "L L"    # ben(e) es-
    foot3 = "L L"    # -se sō-
    foot4 = "L L"    # -lī quom
    foot5 = "S S L"  # si-bi sit
    foot6 = "L S"    # mā-le (final syllable anceps, treated as long for the foot's value)

    # Print the full scansion with feet separated by a space.
    print(f"{foot1} {foot2} {foot3} {foot4} {foot5} {foot6}")

scan_latin_line()