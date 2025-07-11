def solve_language_riddle():
    """
    This function identifies two languages based on a set of orthographic clues
    and prints the result.
    """

    # --- Analysis for Language 'a' ---
    # Clue 1: No 'k' or 'w' in the orthography.
    # Clue 2: Contains the letter 'à'.
    # Clue 3: Is a living, officially recognized language.
    # Deduction: The traditional Italian alphabet has 21 letters, excluding j, k, w, x, and y.
    # These letters are only used for loanwords. The letter 'à' is a standard vowel with a
    # grave accent (e.g., in 'città' or 'perché'). Italian is the official language of Italy.
    # Therefore, language 'a' is identified as Italian.
    language_a = "Italian"

    # --- Analysis for Language 'b' ---
    # Clue 1: Uses the letter combination "ggj" very widely.
    # Clue 2: Uses the letter combination "skt" very widely.
    # Clue 3: Is a living, officially recognized language.
    # Deduction: The trigraph 'ġġ' (often written 'gg' in contexts without the dot)
    # represents a geminated (doubled) sound in Maltese. When followed by 'j' it forms
    # combinations like in 'oġġezzjoni'. Maltese orthography is unique in this regard.
    # The cluster 'skt' is also found in Maltese words like 'jiskrivi' (he writes).
    # Maltese is the official language of Malta.
    # Therefore, language 'b' is identified as Maltese.
    language_b = "Maltese"

    # Print the final answer
    print("Based on the linguistic clues:")
    print(f"Language a is: {language_a}")
    print(f"Language b is: {language_b}")

# Execute the function to display the answer
solve_language_riddle()