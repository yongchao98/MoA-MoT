def analyze_latin_sentence():
    """
    This script analyzes the provided Latin sentence to identify adjectives
    modifying the noun "mariti" based on grammatical rules.
    """
    print("Analyzing which adjectives modify the noun 'mariti'...")

    # The properties of our target noun "mariti" (husband)
    target_noun = {"word": "mariti", "case": "genitive", "number": "singular", "gender": "masculine"}

    # List of adjectives/participles in the sentence and their determined function
    potential_modifiers = [
        {"word": "suppositum", "modifies": "Maxentium"},
        {"word": "muliebri", "modifies": "arte"},
        {"word": "laborantis", "modifies": "mariti"},
        {"word": "gratissimi", "modifies": "partus"},
        {"word": "coepti", "modifies": "partus"},
    ]

    # Find which of the potential modifiers actually modify "mariti"
    modifying_list = []
    for word_info in potential_modifiers:
        if word_info["modifies"] == target_noun["word"]:
            modifying_list.append(word_info["word"])
            print(f"- Found: '{word_info['word']}' modifies 'mariti'.")
        else:
            print(f"- Checked: '{word_info['word']}' modifies '{word_info['modifies']}', not 'mariti'.")


    count = len(modifying_list)

    print("\n--- Final Result ---")
    # This fulfills the prompt's requirement to output numbers in an "equation" format
    print(f"Number of adjectives modifying 'mariti' = {count}")
    if count > 0:
        print(f"The list of modifying adjective(s) is: {modifying_list}")
    else:
        print("The list of modifying adjective(s) is: None")

analyze_latin_sentence()