def analyze_trisyllabic_laxing():
    """
    Analyzes a list of words to determine which has not undergone
    trisyllabic laxing.
    """
    word_analysis = {
        "southern": {
            "base_word": "south",
            "base_vowel": "/aʊ/",
            "derived_vowel": "/ʌ/",
            "syllables": 2,
            "analysis": "Vowel laxing occurred (from /aʊ/ to /ʌ/), but in a two-syllable word. This is related to but not strictly 'trisyllabic' laxing."
        },
        "derivative": {
            "base_word": "derive",
            "base_vowel": "/aɪ/",
            "derived_vowel": "/ɪ/",
            "syllables": 4,
            "analysis": "Classic example of trisyllabic laxing. The long /aɪ/ in 'derive' becomes short /ɪ/ in the antepenultimate syllable of 'derivative'."
        },
        "serenity": {
            "base_word": "serene",
            "base_vowel": "/iː/",
            "derived_vowel": "/ɛ/",
            "syllables": 4,
            "analysis": "Classic example of trisyllabic laxing. The long /iː/ in 'serene' becomes short /ɛ/ in the antepenultimate syllable of 'serenity'."
        },
        "pleasant": {
            "base_word": "please",
            "base_vowel": "/iː/",
            "derived_vowel": "/ɛ/",
            "syllables": 2,
            "analysis": "Vowel laxing occurred (from /iː/ to /ɛ/), but in a two-syllable word. This is pre-cluster shortening, not trisyllabic laxing."
        },
        "gratitude": {
            "base_word": "grateful",
            "base_vowel": "/eɪ/",
            "derived_vowel": "/æ/",
            "syllables": 3,
            "analysis": "Classic example of trisyllabic laxing. The tense vowel /eɪ/ in the root is laxed to /æ/ in the antepenultimate syllable of 'gratitude'."
        },
        "shadow": {
            "base_word": "N/A (not derived from a simpler English base via this process)",
            "base_vowel": "N/A",
            "derived_vowel": "/æ/",
            "syllables": 2,
            "analysis": "Has not undergone this process. Its /æ/ vowel is not the result of shortening a tense vowel from a base like 'shade'. Its phonology developed differently."
        }
    }

    final_answer = None

    print("Analysis of words for Trisyllabic Laxing:")
    print("-" * 50)
    for word, data in word_analysis.items():
        print(f"Word: '{word}'")
        print(f"  - Base Word: {data['base_word']}")
        print(f"  - Vowel Change: {data['base_vowel']} -> {data['derived_vowel']}")
        print(f"  - Analysis: {data['analysis']}")
        if "Has not undergone" in data['analysis']:
            final_answer = word
        print("-" * 50)

    print(f"\nConclusion: All words except '{final_answer}' show a historical vowel shortening related to a base form.")
    print(f"'{final_answer}' is the correct choice because its historical development does not include this specific sound change.")

analyze_trisyllabic_laxing()