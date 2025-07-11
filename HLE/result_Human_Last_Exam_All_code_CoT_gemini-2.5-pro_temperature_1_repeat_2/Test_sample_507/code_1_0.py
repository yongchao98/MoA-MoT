def find_similar_words():
    """
    This function identifies and prints a comparison of words from two
    languages in Asia that have no direct contact but share similar words
    for "mom" and "dad".
    """
    languages = {
        "Korean": {
            "culture": "East Asian",
            "mom": "eomma (엄마)",
            "dad": "appa (아빠)",
            "broom": "bitjaru (빗자루)"
        },
        "Tamil": {
            "culture": "South Asian (Dravidian)",
            "mom": "amma (அம்மா)",
            "dad": "appa (அப்பா)",
            "broom": "vilakkumāṟu (விளக்குமாறு)"
        }
    }

    lang1_name = "Korean"
    lang2_name = "Tamil"
    lang1 = languages[lang1_name]
    lang2 = languages[lang2_name]

    print(f"The two languages are {lang1_name} and {lang2_name}.")
    print(f"These languages are from two different cultures ({lang1['culture']} and {lang2['culture']})")
    print("and have no significant history of direct contact.")
    print("\nHere is a comparison of the words:\n")

    words_to_compare = ["mom", "dad", "broom"]
    for word in words_to_compare:
        print(f"Word: '{word.capitalize()}'")
        print(f"  - {lang1_name}: {lang1[word]}")
        print(f"  - {lang2_name}: {lang2[word]}")
        print("-" * 20)
    
    print("\nNote: The words for 'mom' and 'dad' are nearly identical.")
    print("The word for 'broom' is not similar, but the striking similarity")
    print("of the other two words in these historically non-contacted cultures is a well-known linguistic curiosity.")

find_similar_words()