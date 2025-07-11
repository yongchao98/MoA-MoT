def find_hypothetical_word():
    """
    This function simulates the linguistic evolution of the Old English word for "sister"
    to what it might have become in Modern English without Norse influence.
    """
    # 1. Start with the original Old English word.
    old_english_word = "sweostor"
    print(f"The native Old English word for sister was: {old_english_word}")

    # 2. Apply the first major sound change: the 'sweo' cluster simplifies to 'su'.
    # This was a common path, resulting in the attested Middle English form "suster".
    # old_english_word[0] + old_english_word[3:] gives us "sstor", not what we want.
    # old_english_word.replace("weo", "u") gives us "sustor".
    middle_english_form = "sustor"
    print(f"Step 1: The 'sweo' vowel cluster simplified to 'su'. Transformation: {old_english_word} -> {middle_english_form}")

    # 3. Apply the second major sound change: the unstressed ending '-or' reduces to '-er'.
    # This is a consistent change seen in words like brother (from brōþor).
    final_form = middle_english_form.replace("or", "er")
    print(f"Step 2: The unstressed ending '-or' reduced to '-er'. Transformation: {middle_english_form} -> {final_form}")

    # 4. Present the final result.
    print(f"\nConclusion: Without the influence from Old Norse 'systir', the native word '{old_english_word}' would have most likely evolved into '{final_form}'.")


find_hypothetical_word()