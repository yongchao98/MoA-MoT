def find_native_english_word():
    """
    This function traces the hypothetical evolution of the Old English word for "sister"
    to determine what its Modern English form would be without Norse influence.
    """

    # Step 1: Identify the current word and its non-native origin.
    current_word = "sister"
    norse_origin = "systir"
    print(f"Analyzing the word '{current_word}':")
    print(f"- The current word was borrowed directly from the Old Norse '{norse_origin}' during the Viking Age.")
    print("- To find the native word, we must ignore this borrowing and look back to Old English.")
    print("\n" + "="*50 + "\n")

    # Step 2: Identify the original Old English (Anglo-Saxon) word.
    old_english_word = "sweostor"
    print(f"Finding the native root:")
    print(f"- The original Old English word for sister was '{old_english_word}'.")
    print("- This is a cognate of German 'Schwester' and Dutch 'zuster'.")
    print("\n" + "="*50 + "\n")

    # Step 3: Trace the likely linguistic evolution to Modern English.
    middle_english_form = "suster"
    predicted_word = "suster"
    print("Projecting the evolution without Norse influence:")
    print(f"- Old English '{old_english_word}' would have naturally simplified.")
    print(f"- By the Middle English period, a very common form was '{middle_english_form}'.")
    print(f"- The short 'u' sound in '{middle_english_form}' would follow the same path as the 'u' in 'sunu' (son) or 'lufu' (love).")
    print("- This vowel sound would become the modern /ʌ/ sound (the 'uh' in 'duster').")
    print("\n" + "="*50 + "\n")

    # Step 4: State the final conclusion.
    print("Conclusion:")
    print(f"Therefore, had the Norse never invaded, the Modern English word for '{current_word}'")
    print("would most likely be:")
    print(f"\n>>> {predicted_word}")


# Run the analysis
find_native_english_word()