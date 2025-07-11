def find_critical_word_position():
    """
    Analyzes a passage to find the word position where a metonymic interpretation
    would cause elevated reading times, heading off a garden path.
    """
    passage = "Experience will be a decisive factor in this battle. The ships will outmaneuver the enemy fleet because the old man the boats."

    # Split the passage into words. Punctuation is kept with the word as in the example.
    words = passage.split()

    # The critical word is "old", where the metonymic interpretation ("the old" -> "the old people") occurs.
    critical_word = "old"

    try:
        # Find the index of the critical word. list.index() is 0-based.
        word_index = words.index(critical_word)
        # Ordinal position is 1-based.
        ordinal_position = word_index + 1

        print(f"The full passage is: \"{passage}\"")
        print(f"The list of words (tokens) is: {words}")
        print("-" * 20)
        print(f"The critical word where the metonymic interpretation occurs is '{critical_word}'.")
        # As per the instructions, output the number in the final result.
        print(f"The final equation is: Ordinal Position = Index + 1 => {ordinal_position} = {word_index} + 1")
        print(f"The word '{critical_word}' is at ordinal position: {ordinal_position}")

    except ValueError:
        print(f"The critical word '{critical_word}' was not found in the passage.")

if __name__ == "__main__":
    find_critical_word_position()