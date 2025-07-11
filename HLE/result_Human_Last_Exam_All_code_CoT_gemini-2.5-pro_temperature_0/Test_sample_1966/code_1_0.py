def solve_shakespeare_dante_question():
    """
    Determines which Shakespearean title characters are mentioned in The Divine Comedy.
    """
    # All unique characters from the answer choices
    all_characters = {
        "Julius Caesar", "Pericles", "Cleopatra", "King John", "Troilus", "Antony"
    }

    # Characters confirmed to be mentioned by name in The Divine Comedy
    # - Julius Caesar: Inferno, Canto IV
    # - Cleopatra: Inferno, Canto V
    dante_mentions = {"Julius Caesar", "Cleopatra"}

    print("Checking which characters from the list appear in The Divine Comedy...")
    print("-" * 60)

    found_characters = []
    for char in sorted(list(all_characters)):
        if char in dante_mentions:
            print(f"Result for '{char}': Found")
            found_characters.append(char)
        else:
            print(f"Result for '{char}': Not Found")

    print("-" * 60)
    print("The Shakespearean title characters mentioned by name in The Divine Comedy are:")
    # Sort for consistent ordering before printing
    final_list = sorted(found_characters)
    print(", ".join(final_list))

    # Match the result with the given answer choices
    answer_choices = {
        "A": ["Julius Caesar", "Pericles"],
        "B": ["Julius Caesar", "Cleopatra", "King John"],
        "C": ["Julius Caesar", "Troilus", "Antony"],
        "D": ["Julius Caesar", "Cleopatra"],
        "E": ["Julius Caesar", "Antony", "Cleopatra"]
    }

    correct_choice = "Unknown"
    for choice, names in answer_choices.items():
        # Compare sets for an order-independent match
        if set(names) == set(final_list):
            correct_choice = choice
            break

    print(f"\nThis list matches answer choice: {correct_choice}")


solve_shakespeare_dante_question()