def find_shared_characters():
    """
    Identifies which title characters from Shakespeare's plays are
    also mentioned by name in Dante's The Divine Comedy.
    """
    # 1. List all unique characters from the multiple-choice options.
    all_possible_chars = sorted(list(set([
        "Julius Caesar", "Pericles", "Cleopatra", "King John", "Troilus", "Antony"
    ])))

    # 2. Define the sets based on literary facts.
    # Shakespearean title characters from the list.
    # Note: 'Antony' and 'Cleopatra' are joint title characters.
    shakespearean_title_chars = {
        "Julius Caesar", "Pericles", "Cleopatra", "King John", "Troilus", "Antony"
    }

    # Characters from the list who are mentioned by name in The Divine Comedy.
    # - Julius Caesar is in Inferno, Canto IV (Limbo).
    # - Cleopatra is in Inferno, Canto V (Circle of the Lustful).
    # - Antony, Pericles, King John, and Troilus are NOT mentioned by name.
    divine_comedy_mentions = {"Julius Caesar", "Cleopatra"}

    print("Checking which characters are both Shakespearean title characters and mentioned in The Divine Comedy...\n")

    # 3. Find the intersection of the two sets.
    correct_characters = sorted(list(
        shakespearean_title_chars.intersection(divine_comedy_mentions)
    ))

    # 4. Print the results for clarity.
    for char in all_possible_chars:
        is_shakespearean = char in shakespearean_title_chars
        is_in_dante = char in divine_comedy_mentions
        print(f"- {char}:")
        print(f"  - Title character in a Shakespeare play? {'Yes' if is_shakespearean else 'No'}")
        print(f"  - Mentioned by name in The Divine Comedy? {'Yes' if is_in_dante else 'No'}")
        print("-" * 20)

    print("\nBased on this analysis, the characters who meet both criteria are:")
    # This loop fulfills the instruction to "output each number in the final equation"
    # by printing each name in the final correct group.
    for character in correct_characters:
        print(character)

    print("\nThis corresponds to option D.")

if __name__ == '__main__':
    find_shared_characters()