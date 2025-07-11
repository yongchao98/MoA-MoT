import collections

def solve_shakespeare_dante_question():
    """
    This function identifies Shakespearean title characters who are also mentioned
    by name in Dante Alighieri's "The Divine Comedy".
    """
    # Step 1: Define the set of unique Shakespearean title characters from the answer choices.
    shakespeare_title_characters_from_options = {
        "Julius Caesar",
        "Pericles",
        "Cleopatra",
        "King John",
        "Troilus",
        "Antony"
    }

    # Step 2: Define the set of characters from the above list who are explicitly
    # named in "The Divine Comedy", based on literary analysis.
    characters_named_in_divine_comedy = {
        "Julius Caesar",  # Found in Inferno, Canto IV
        "Cleopatra"       # Found in Inferno, Canto V
    }

    # Step 3: Find the intersection of the two sets to identify the common characters.
    common_characters = shakespeare_title_characters_from_options.intersection(characters_named_in_divine_comedy)

    # Step 4: Sort the list for consistent output and print the final result.
    final_list = sorted(list(common_characters))

    print("The following Shakespearean title characters are mentioned by name in The Divine Comedy:")
    # As requested, printing each name in the final list.
    print(", ".join(final_list))

solve_shakespeare_dante_question()
<<<D>>>