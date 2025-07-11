def find_common_characters():
    """
    Finds which Shakespearean title characters from the options
    are also mentioned by name in Dante's 'The Divine Comedy'.
    """

    # Title characters from Shakespeare plays mentioned in the answer choices.
    shakespeare_title_characters = {
        "Julius Caesar",
        "Pericles",
        "Cleopatra",
        "King John",
        "Troilus",
        "Antony"
    }

    # A set of characters confirmed to be mentioned by name in 'The Divine Comedy'.
    # Note: While a statesman named Pericles is in 'The Divine Comedy', he is not the
    # title character of Shakespeare's play 'Pericles, Prince of Tyre'.
    # Mark Antony is not mentioned by name.
    divine_comedy_characters = {
        "Julius Caesar",
        "Cleopatra",
        # Adding some other examples for context
        "Hector",
        "Aeneas",
        "Homer",
        "Aristotle",
        "Plato",
        "Brutus",
        "Cassius",
        "Helen of Troy",
        "Paris",
        "Tristan"
    }

    # Find the intersection of the two sets
    common_characters = shakespeare_title_characters.intersection(divine_comedy_characters)

    print("The Shakespearean title characters mentioned by name in 'The Divine Comedy' are:")
    for char in sorted(list(common_characters)):
        print(f"- {char}")

find_common_characters()