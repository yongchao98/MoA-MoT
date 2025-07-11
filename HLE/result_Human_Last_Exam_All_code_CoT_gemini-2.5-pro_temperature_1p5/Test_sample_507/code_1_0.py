import textwrap

def find_similar_language_words():
    """
    Identifies and explains two Asian languages from different cultures
    with no direct contact that share surprisingly similar words for
    "broom," "mom," and "dad".
    """

    lang1 = "Turkish"
    lang2 = "Malay"
    culture1 = "Turkey (West Asia)"
    culture2 = "Malaysia/Indonesia (Southeast Asia)"

    explanation = (
        f"A compelling example of coincidental linguistic similarities can be found between {lang1} and {lang2}. These languages belong to unrelated families (Turkic and Austronesian, respectively) and their native cultures developed thousands of miles apart in {culture1} and {culture2} with no significant direct historical contact."
    )

    # Word data
    word_pairs = {
        "Dad": {
            "lang1_word": "baba",
            "lang2_word": "bapa",
            "comment": "These words are nearly identical in sound and meaning."
        },
        "Mom": {
            "lang1_word": "ana",
            "lang2_word": "emak",
            "comment": "While not identical, both are simple words for 'mother' that prominently feature an 'a' vowel and a nasal consonant ('n'/'m')."
        },
        "Broom": {
            "lang1_word": "süpürge (from the verb süpürmek, 'to sweep')",
            "lang2_word": "penyapu (from the root word sapu, 'to sweep')",
            "comment": "The root words for the act of sweeping, 'süpür-' and 'sapu', are phonetically similar, both starting with 's' and featuring a 'p' sound. This is a striking coincidence for a non-nursery word."
        }
    }

    print("--- Linguistic Puzzle Answer ---")
    print(textwrap.fill(explanation, 80))
    print("\n" + "="*40 + "\n")
    
    print("Here is the word-by-word comparison:\n")

    for term, data in word_pairs.items():
        print(f"WORD: \"{term}\"")
        print(f"  - In {lang1}: {data['lang1_word']}")
        print(f"  - In {lang2}: {data['lang2_word']}")
        print(f"  - Similarity: {data['comment']}\n")

# Execute the function to display the answer
find_similar_language_words()