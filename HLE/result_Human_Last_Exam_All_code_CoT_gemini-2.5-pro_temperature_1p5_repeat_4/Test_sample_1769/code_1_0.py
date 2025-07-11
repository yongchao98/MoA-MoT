def find_hypothetical_word():
    """
    This function explains and calculates the hypothetical Modern English word for 'sister'
    if the Norse influence had never occurred.
    """
    print("Analyzing the etymology to find the hypothetical word for 'sister'...")
    print("-" * 60)

    # Step 1: Explain the origin of the current word
    current_word = "sister"
    old_norse_origin = "systir"
    print(f"The current word, '{current_word}', is not native to Old English.")
    print(f"It was borrowed from the Old Norse word '{old_norse_origin}' during the Viking Age (roughly 8th to 11th centuries).\n")

    # Step 2: Identify the original Old English word
    old_english_word = "sweostor"
    print(f"The original Old English (Anglo-Saxon) word for sister was '{old_english_word}'.\n")

    # Step 3: Trace the linguistic evolution
    print("Without the Norse influence, '{old_english_word}' would have evolved naturally.")
    print("Similar to how other Old English words changed:")
    print(" - 'mōdor' became 'mother'")
    print(" - 'brōþor' became 'brother'")
    print(f"Following these patterns, the Old English '{old_english_word}' would have likely simplified.\n")

    # Step 4: Present the final predicted word
    predicted_word = "swester"
    print(f"The word '{old_english_word}' would have evolved into the Middle English 'swester', and likely would have remained as such.")
    print("In fact, the word 'swester' existed in English for a time but was eventually outcompeted by the Norse version.")
    print("-" * 60)
    print("Therefore, the hypothetical Modern English word would be:")
    print(f"\n>>> {predicted_word} <<<")


if __name__ == "__main__":
    find_hypothetical_word()