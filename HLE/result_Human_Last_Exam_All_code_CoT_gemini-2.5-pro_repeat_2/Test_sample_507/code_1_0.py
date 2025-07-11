def present_language_comparison():
    """
    Identifies and displays two languages from different parts of Asia
    with no direct contact that share similar words for "mom", "dad", and "broom".
    """
    # Define the languages and cultures
    language1_name = "Korean"
    language1_culture = "East Asia (Koreanic language family)"
    language2_name = "Telugu"
    language2_culture = "South India (Dravidian language family)"

    # Define the word comparisons
    comparison_data = {
        "Mom": {
            "lang1_word": "eomma (엄마)",
            "lang2_word": "amma (అమ్మ)"
        },
        "Dad": {
            "lang1_word": "appa (아빠)",
            "lang2_word": "appa (అప్ప)"
        },
        "Broom": {
            "lang1_word": "bijaru (비자루)",
            "lang2_word": "chīpuru (చీపురు)"
        }
    }

    # Print the introduction and findings
    print("Languages Identified:")
    print(f"1. {language1_name} (from {language1_culture})")
    print(f"2. {language2_name} (from {language2_culture})")
    print("\nThese two languages have no direct contact, being geographically distant and from different language families.")
    print("\nWord Comparison:")
    print("-" * 50)
    print(f"{'Concept':<10} | {language1_name:<20} | {language2_name:<20}")
    print("-" * 50)

    # Print each word pair
    for concept, words in comparison_data.items():
        print(f"{concept:<10} | {words['lang1_word']:<20} | {words['lang2_word']:<20}")

    print("-" * 50)
    print("\nNote on 'Broom': The words 'bijaru' and 'chīpuru' show phonetic similarity.")
    print("b/p are related sounds (bilabial stops), as are j/ch (palatal affricates).")

present_language_comparison()