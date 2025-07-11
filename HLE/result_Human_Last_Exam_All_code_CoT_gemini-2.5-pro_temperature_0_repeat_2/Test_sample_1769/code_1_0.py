def find_hypothetical_word():
    """
    This function determines what the Modern English word for "sister" would be
    without the influence of the Old Norse language from Viking invasions.
    """

    # 1. The current word "sister" comes from Middle English "sister", which was
    #    heavily influenced by, and largely replaced, the native word due to the
    #    Old Norse cognate "systir".
    old_norse_word = "systir"
    current_word = "sister"

    # 2. The original, native Old English (Anglo-Saxon) word was "sweostor".
    old_english_word = "sweostor"

    # 3. To find the modern equivalent, we trace the likely evolution of "sweostor"
    #    through Middle English to Modern English.

    #    - From Old English to Middle English: The 'eo' sound in 'sweostor'
    #      would typically simplify to a short 'e' sound. The unstressed 'o' in the
    #      second syllable would also reduce to an 'e'.
    #      This gives us the Middle English form "swester". This form actually
    #      existed and was used alongside "sister".
    middle_english_word = "swester"

    #    - From Middle English to Modern English: The short 'e' vowel in "swester"
    #      would have been stable and not undergone the Great Vowel Shift in the
    #      same way a long vowel would have. Its pronunciation would remain,
    #      similar to the evolution of "fether" to "feather".
    hypothetical_modern_word = "swester"

    print(f"The original Old English word was: '{old_english_word}'")
    print(f"It was replaced by the Old Norse-influenced word: '{current_word}' (from O.N. '{old_norse_word}')")
    print("-" * 30)
    print("Without Norse influence, the native word would have evolved as follows:")
    print(f"Old English '{old_english_word}' -> Middle English '{middle_english_word}' -> Modern English '{hypothetical_modern_word}'")
    print("-" * 30)
    print(f"Therefore, the Modern English word for 'sister' would likely be: {hypothetical_modern_word}")

find_hypothetical_word()