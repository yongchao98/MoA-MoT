def find_similar_words():
    """
    This function identifies and prints information about two languages from
    different Asian cultures with no direct contact that share similar words for
    "mom", "dad", and "broom".
    """

    # Data structure holding the information
    language_data = {
        "Korean": {
            "culture": "East Asia (Korean Peninsula)",
            "words": {
                "Mom": "Eomma (엄마)",
                "Dad": "Appa (아빠)",
                "Broom": "Bitjaru (빗자루)"
            }
        },
        "Telugu": {
            "culture": "South India (Andhra Pradesh & Telangana states)",
            "words": {
                "Mom": "Amma (అమ్మ)",
                "Dad": "Appa (అప్ప)",
                "Broom": "Cheepuru (చీపురు)"
            }
        }
    }

    lang1_name = "Korean"
    lang2_name = "Telugu"
    lang1 = language_data[lang1_name]
    lang2 = language_data[lang2_name]

    # --- Introduction ---
    print(f"The two languages are {lang1_name} and {lang2_name}.\n")
    print(f"{lang1_name}: A Koreanic language spoken in {lang1['culture']}.")
    print(f"{lang2_name}: A Dravidian language spoken in {lang2['culture']}.")
    print("\nThese cultures have no significant direct historical contact, making the similarities a remarkable coincidence.")
    print("-" * 60)

    # --- Word Comparison ---
    words_to_compare = ["Mom", "Dad", "Broom"]
    for word_type in words_to_compare:
        print(f"\nComparing the word for '{word_type}':")
        print(f"  In {lang1_name}, it is: {lang1['words'][word_type]}")
        print(f"  In {lang2_name}, it is: {lang2['words'][word_type]}")

# Execute the function to print the result
find_similar_words()
