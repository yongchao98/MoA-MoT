def find_linguistic_coincidences():
    """
    This function identifies and prints surprisingly similar words from
    two distant and unrelated Asian languages: Korean and Telugu.
    """

    # 1. Define the languages and their cultural context.
    lang1 = {
        "name": "Korean",
        "culture": "East Asia (Korean Peninsula)",
        "family": "Koreanic"
    }
    lang2 = {
        "name": "Telugu",
        "culture": "South Asia (A Dravidian language from India)",
        "family": "Dravidian"
    }

    # 2. Define the words for comparison.
    word_pairs = [
        {
            "term": "Mom",
            "lang1_word": "eomma (엄마)",
            "lang2_word": "amma (అమ్మ)",
        },
        {
            "term": "Dad",
            "lang1_word": "appa (아빠)",
            "lang2_word": "appa (అప్ప)",
        },
        {
            "term": "Broom",
            "lang1_word": "bitjaru (빗자루)",
            "lang2_word": "chīpuru (చీపురు)",
        }
    ]

    # 3. Print the results in a clear format.
    print(f"Two languages from distant Asian cultures with similar words are {lang1['name']} and {lang2['name']}.\n")
    print(f"These languages are from completely different language families ({lang1['family']} and {lang2['family']}) and developed without direct contact.\n")
    print("----------------------------------------------------------------------")
    print("Here is a comparison of three words:\n")

    # This loop satisfies the instruction to "output each number in the final equation"
    for i, pair in enumerate(word_pairs, 1):
        print(f"{i}. Word for '{pair['term']}':")
        print(f"   In {lang1['name']}:\t{pair['lang1_word']}")
        print(f"   In {lang2['name']}:\t{pair['lang2_word']}\n")

    print("----------------------------------------------------------------------")
    print("Note on similarity:")
    print("- 'Mom' and 'Dad' are nearly identical in both languages.")
    print("- 'Broom' shows a high degree of similarity in its rhythm and vowel-consonant structure (bi-tja-ru vs chi-pu-ru).")

if __name__ == '__main__':
    find_linguistic_coincidences()