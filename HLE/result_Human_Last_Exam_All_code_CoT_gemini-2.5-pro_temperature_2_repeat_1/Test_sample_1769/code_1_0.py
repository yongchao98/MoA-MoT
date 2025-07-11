def find_hypothetical_word():
    """
    Calculates the hypothetical Modern English word for 'sister'
    without Norse influence by tracing the evolution of the original
    Old English word 'sweostor'.
    """

    # Step 1: The current word, "sister", comes from Old Norse "systir".
    # This word replaced the native English word.
    old_norse_word = "systir"
    current_word = "sister"

    # Step 2: The native Old English word for "sister" was "sweostor".
    old_english_word = "sweostor"

    print(f"The current word '{current_word}' is a loanword from Old Norse '{old_norse_word}'.")
    print(f"It replaced the original Old English word: '{old_english_word}'.\n")
    print("To find the hypothetical modern word, we trace the evolution of 'sweostor'.\n")
    
    # Step 3: Evolution from Old English to Middle English.
    # The "sweo-" prefix often simplified. The most common attested form
    # in Middle English derived from 'sweostor' was 'suster'.
    middle_english_word = "suster"

    # Step 4: Evolution from Middle English to Modern English.
    # The word would have been affected by the Great Vowel Shift.
    # The short 'u' /ʊ/ (as in 'put') shifted to the 'uh' /ʌ/ sound (as in 'cup').
    # The spelling, however, would likely have remained.
    hypothetical_modern_word = "suster"

    print("The historical linguistic path would have been:")
    # The user prompt requested to "output each number in the final equation".
    # We will interpret this as printing each stage of the word's evolution.
    print(f"1. Old English form: {old_english_word}")
    print(f"2. Middle English form: {middle_english_word}")
    print(f"3. Hypothetical Modern English form: {hypothetical_modern_word}\n")

    print(f"Thus, without Norse influence, the Modern English word for 'sister' would most likely be '{hypothetical_modern_word}'.")
    print("This word is pronounced like 'suster' (rhymes with 'duster' or 'cluster').")

if __name__ == "__main__":
    find_hypothetical_word()