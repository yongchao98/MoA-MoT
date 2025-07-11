def find_shared_characters():
    """
    Finds which Shakespearean title characters from a given list
    are also mentioned by name in Dante's "The Divine Comedy".
    """
    # Step 1: List of Shakespearean title characters from the answer choices.
    shakespeare_chars = sorted(list(set(['Julius Caesar', 'Pericles', 'Cleopatra', 'King John', 'Troilus', 'Antony'])))

    # Step 2: A list of figures explicitly mentioned by name in "The Divine Comedy".
    # This is based on literary analysis of Dante's text.
    # - Julius Caesar is in Inferno, Canto IV (Limbo).
    # - Cleopatra is in Inferno, Canto V (Circle of the Lustful).
    # - Troilus is generally believed to be "the good Trojan" ("il buon Troyano")
    #   in Inferno, Canto IV, but is not mentioned by his proper name.
    # - Antony (Mark Antony), Pericles, and King John are not mentioned by name.
    divine_comedy_chars = ['Julius Caesar', 'Cleopatra', 'Hector', 'Aeneas', 'Homer', 'Socrates', 'Plato', 'Saladin', 'Helen', 'Paris', 'Tristan']

    print("Checking which Shakespearean title characters are in The Divine Comedy...")
    print("-" * 60)

    found_characters = []
    for char in shakespeare_chars:
        if char in divine_comedy_chars:
            found_characters.append(char)
            location = ""
            if char == "Julius Caesar":
                location = "Inferno, Canto IV (Limbo)"
            elif char == "Cleopatra":
                location = "Inferno, Canto V (Circle of the Lustful)"
            print(f"- Found: {char} (mentioned in {location})")
        else:
            print(f"- Not Found: {char}")

    print("-" * 60)
    print("Based on this analysis, the title characters from Shakespeare's plays who are also mentioned by name in The Divine Comedy are:")
    # The final print uses a formatted string to show the result clearly.
    final_result = ", ".join(found_characters)
    print(final_result)


find_shared_characters()