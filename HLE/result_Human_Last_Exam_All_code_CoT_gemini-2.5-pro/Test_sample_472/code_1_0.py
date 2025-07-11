def scan_latin_line():
    """
    Scans the provided Latin line and prints the metrical scansion.

    The line is: "et tibi bene esse soli quom sibi sit male"
    The meter is a modified trochaic septenarius.

    Scansion breakdown:
    1. et tibi: Scanned as a dactyl (L S S). 'et' is long by position.
    2. bene esse: Elides to 'benesse'. Scanned as S L S. Due to iambic shortening,
       the initial S L ('be-nes') becomes S S. The foot becomes a tribrach (S S S).
    3. soli: Scanned as a trochee (L S). 'o' is long by nature.
    4. quom: A single long syllable, often marking the caesura (mid-line pause). (L)
    5. sibi sit: Scanned as S S L. Due to iambic shortening, the S L ('bi sit')
       becomes S S. The foot becomes a tribrach (S S S).
    6. male: Scanned as an iamb (S L), with the final syllable long by convention at the end of a verse.
    """

    # The final scansion is derived from the steps above.
    # Feet are: (L S S) (S S S) (L S) (L) (S S S) (S L)
    foot1 = "L S S"
    foot2 = "S S S"
    foot3 = "L S"
    foot4 = "L"
    foot5 = "S S S"
    foot6 = "S L"

    final_scansion = f"{foot1} {foot2} {foot3} {foot4} {foot5} {foot6}"
    print(final_scansion)

scan_latin_line()