def analyze_word_avalanche():
    """
    Analyzes choices for a "True Avalanche" based on a set of criteria.
    A True Avalanche is a pun where syllables are repeated phonetically.
    """

    # The required word and description for the pun
    target_word = "computer"
    syllable_pattern = "com-pu-ter"
    description = "My software tells the birds when and where to relieve themselves."

    # The answer choices provided
    choices = {
        'A': "Computers comp few urns",
        'B': "Computer: \"Come poo, tern!\"",
        'C': "Computer: \"Come, pee you turkey!\"",
        'D': "Comb pewter computer",
        'E': "Computer needs bird relieves"
    }

    print(f"Analyzing choices for a True Avalanche for the word '{target_word}'...\n")
    print(f"The target syllable pattern is: {syllable_pattern}")
    print(f"The description to match is: '{description}'\n")

    # Store analysis results
    analysis = {
        'A': {
            "pun": "comp few urns",
            "pun_quality": "Weak phonetic match.",
            "description_match": "No. Does not relate to telling a bird anything."
        },
        'B': {
            "pun": "\"Come poo, tern!\"",
            "pun_quality": "Excellent phonetic match to 'com-pu-ter'.",
            "description_match": "Yes. It's a command to a 'tern' (a bird) to 'poo'."
        },
        'C': {
            "pun": "\"Come, pee you turkey!\"",
            "pun_quality": "Weak phonetic match. Only 'Come' matches 'com'.",
            "description_match": "Yes, it matches the theme, but the pun is poor."
        },
        'D': {
            "pun": "Comb pewter",
            "pun_quality": "Excellent phonetic match to 'computer'.",
            "description_match": "No. 'Comb pewter' is unrelated to the description."
        },
        'E': {
            "pun": "N/A",
            "pun_quality": "Not a pun or word avalanche.",
            "description_match": "No, this is a statement about the theme, not an example of the pun."
        }
    }

    best_choice = None
    for key, text in choices.items():
        if analysis[key]["pun_quality"].startswith("Excellent") and analysis[key]["description_match"].startswith("Yes"):
            best_choice = key
            break

    print("--- Analysis Results ---")
    for key, details in analysis.items():
        print(f"Choice {key}: {choices[key]}")
        print(f"  - Pun Quality: {details['pun_quality']}")
        print(f"  - Description Match: {details['description_match']}\n")


    if best_choice:
        print("--- Conclusion ---")
        print(f"The best choice is B.")
        print("It is the only option that is both an excellent phonetic pun AND perfectly matches the provided description.")
        print(f"\nThe final answer is: {choices[best_choice]}")
    else:
        print("No choice perfectly fits all criteria.")

# Run the analysis
analyze_word_avalanche()