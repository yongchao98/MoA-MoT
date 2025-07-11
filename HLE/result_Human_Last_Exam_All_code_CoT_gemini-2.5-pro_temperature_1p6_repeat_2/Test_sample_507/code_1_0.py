def find_similar_words():
    """
    Identifies and prints two languages from different Asian cultures with
    no direct contact that have similar words for 'broom', 'mom', and 'dad'.
    """

    language1 = "Tagalog (Philippines)"
    language2 = "Mongolian (Mongolia)"

    word_pairs = {
        "Dad": {
            "lang1_word": "Ama",
            "lang2_word": "Aav (Аав)",
            "similarity_note": "Both start with 'A' and a labial consonant (m/v)."
        },
        "Mom": {
            "lang1_word": "Ina",
            "lang2_word": "Eej (Ээж)",
            "similarity_note": "Both are vowel-initial words with a similar feel, common in 'nursery words'."
        },
        "Broom": {
            "lang1_word": "Walis",
            "lang2_word": "Bilis (билис)",
            "similarity_note": "A surprising phonetic coincidence in structure and sound (w/b, l, s)."
        }
    }

    print("The two languages are Tagalog and Mongolian.\n")
    print("These represent two distinct Asian cultures (Philippines and Mongolia) from different language families (Austronesian and Mongolic) with no direct historical contact.\n")
    print(f"{'Concept':<10} | {language1:<25} | {language2:<25} | {'Note':<60}")
    print("-" * 125)

    for concept, words in word_pairs.items():
        lang1_word = words['lang1_word']
        lang2_word = words['lang2_word']
        note = words['similarity_note']
        print(f"{concept:<10} | {lang1_word:<25} | {lang2_word:<25} | {note}")

find_similar_words()