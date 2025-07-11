def find_scansion_issue():
    """
    Analyzes a sestina to find a word that breaks its structural rules.
    """
    # Step 1: Define the six key end-words from the first stanza.
    # The end-words are: vainly, fly, call, tea, birds, past.
    # We will accept simple variations like 'calls' for 'call' and 'birds' for 'bird'.
    key_word_roots = ['vainly', 'fly', 'call', 'tea', 'bird', 'past']

    # Step 2: Examine the final lines of the poem (the envoi).
    # A sestina's envoi must use the key words. Let's check the end-words here.
    final_lines = [
        "fly like birds past dark high-windowed rooms.",
        "An empty teacup rattles. No-one calls."
    ]

    print("Analyzing the poem's structure to find the issue...")
    print("A sestina is defined by the repetition of six key end-words from the first stanza.")
    print(f"The key words for this poem are: {', '.join(key_word_roots)}.\n")

    # Step 3: Check each final line for a valid end-word.
    problem_found = False
    for line in final_lines:
        # Extract the last word, removing punctuation.
        last_word = line.strip().rstrip('.').split()[-1]

        # Check if the last word contains any of the key word roots.
        is_valid_end_word = any(root in last_word.lower() for root in key_word_roots)

        if not is_valid_end_word:
            problem_found = True
            print(f"Error found in the line: '{line}'")
            print(f"The line ends with the word '{last_word}'.")
            print(f"This word is not one of the six required key words, which breaks the sestina's form.")
            print("\nThe word that causes the scansion (structural) issue is therefore:")
            # Use a print statement to output the final answer as requested
            print(last_word)
            break
    
    if not problem_found:
        print("No structural issues were found in the final lines.")

find_scansion_issue()