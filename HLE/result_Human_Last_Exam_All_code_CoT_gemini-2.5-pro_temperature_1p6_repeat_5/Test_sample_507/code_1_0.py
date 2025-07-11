import collections

def find_and_print_similar_words():
    """
    This script presents two languages from different Asian cultures
    with no direct contact that share surprisingly similar words for
    "mom", "dad", and potentially "broom".
    """

    # Using an ordered dictionary to maintain the presentation order
    language_data = collections.OrderedDict()

    language_data["Korean"] = {
        "Culture/Region": "East Asia",
        "Language Family": "Koreanic",
        "Word for 'mom'": "엄마 (eomma)",
        "Word for 'dad'": "아빠 (appa)",
        "Word for 'broom'": "비 (bi) - This is the root word for broom, found in 빗자루 (bitjaru)."
    }

    language_data["Tamil"] = {
        "Culture/Region": "South Asia (India, Sri Lanka)",
        "Language Family": "Dravidian",
        "Word for 'mom'": "அம்மா (amma)",
        "Word for 'dad'": "அப்பா (appa)",
        "Word for 'broom'": "பிரி (piiri) - This word means 'to separate/comb'. Its phonetic similarity to Korean 'bi' is a subject of linguistic curiosity, although the common word for broom is different."
    }

    print("--- Language Comparison ---\n")

    for lang, data in language_data.items():
        print(f"Language: {lang}")
        print(f"  - Culture/Region: {data['Culture/Region']}")
        print(f"  - Language Family: {data['Language Family']}")
        print("\n  Word Comparisons:")
        # The prompt mentioned "output each number in the final equation!", which seems to be a slight error.
        # Interpreting this as "output each word comparison clearly".
        print(f"  - Mom: {data['Word for \'mom\'']}")
        print(f"  - Dad: {data['Word for \'dad\'']}")
        print(f"  - Broom: {data['Word for \'broom\'']}")
        print("-" * 25)

if __name__ == '__main__':
    find_and_print_similar_words()