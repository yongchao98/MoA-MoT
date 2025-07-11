def find_similar_words():
    """
    Identifies and prints two languages from different Asian cultures
    with similar words for "mom," "dad," and "broom."
    """
    lang1 = "Tagalog (Philippines)"
    lang2 = "Tamil (India, Sri Lanka)"

    print(f"The two languages are {lang1} and {lang2}.")
    print("These languages belong to the Austronesian and Dravidian families, respectively,")
    print("and their cultures developed with no significant direct contact.\n")
    print("Here are the similar words for 'mom', 'dad', and 'broom':\n")

    # Word comparisons
    words = {
        "Concept": ["Tagalog Word", "Tamil Word", "Similarity"],
        "Mom": ["Nanay", "Annai (அன்னை)", "Both are two-syllable, start with 'Na'/'An', and end with 'ay'/'ai' sound."],
        "Dad": ["Tatay", "Thanthai (தந்தை)", "Both are two-syllable, start with 'Ta'/'Tha', and have similar dental stop sounds."],
        "Broom": ["Walis", "Vāriyal (வாரியல்)", "Both start with a W/V sound. 'Walis' is for sweeping, while 'Vāriyal' derives from the verb 'vāri' (to gather/scoop up)."]
    }

    # Print table header
    header = words["Concept"]
    print(f"{header[0]:<15} | {header[1]:<20} | {header[2]:<25}")
    print("-" * 70)

    # Print rows
    for concept in ["Mom", "Dad", "Broom"]:
        row = words[concept]
        print(f"{concept:<15} | {row[0]:<20} | {row[1]:<20}")

    print("\n" + "="*70)
    print("Final words expressing the similarity:")
    print("Nanay / Annai")
    print("Tatay / Thanthai")
    print("Walis / Vāriyal")

find_similar_words()