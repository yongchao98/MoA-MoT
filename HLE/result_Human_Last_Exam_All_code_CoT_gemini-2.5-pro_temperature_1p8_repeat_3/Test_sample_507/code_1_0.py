import collections

def find_and_print_similar_words():
    """
    Identifies and prints a comparison of similar words for "dad", "mom", and "broom"
    in two distant Asian languages: Turkish and Malay.
    """

    # Using an ordered dictionary to maintain the order of words for printing
    word_comparisons = collections.OrderedDict()

    word_comparisons["Dad"] = {
        "language1": "Turkish",
        "word1": "Baba",
        "language2": "Malay",
        "word2": "Bapa"
    }
    word_comparisons["Mom"] = {
        "language1": "Turkish",
        "word1": "Anne",
        "language2": "Malay",
        "word2": "Emak"
    }
    # For "broom", the verb "to sweep" shows the root similarity more clearly.
    # Turkish "süpürge" (broom) comes from the verb "süpürmek" (to sweep).
    # Malay "penyapu" (broom) comes from the verb "sapu" (to sweep).
    word_comparisons["Broom (from the verb 'to sweep')"] = {
        "language1": "Turkish",
        "word1": "süpür- (root of süpürge)",
        "language2": "Malay",
        "word2": "sapu (root of penyapu)"
    }

    print("Two languages from distant Asian cultures (Turkish from West Asia and Malay from Southeast Asia) with similar words:")
    print("="*80)

    for english_word, details in word_comparisons.items():
        print(f"The English word is: '{english_word}'")
        print(f"  - In {details['language1']}, the word is: '{details['word1']}'")
        print(f"  - In {details['language2']}, the word is: '{details['word2']}'")
        print("-" * 40)

# Execute the function to print the comparison.
find_and_print_similar_words()