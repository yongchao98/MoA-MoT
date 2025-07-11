def find_word_without_laxing():
    """
    Analyzes a list of words to find the one that has not undergone
    vowel laxing (shortening) during its historical development.
    """
    word_analysis = [
        {
            "word": "southern",
            "root": "south",
            "root_vowel": "/aʊ/ (tense diphthong)",
            "derived_vowel": "/ʌ/ (lax vowel)",
            "explanation": "The tense diphthong in 'south' becomes a lax vowel in 'southern'. Laxing has occurred.",
            "has_laxed": True
        },
        {
            "word": "derivative",
            "root": "derive",
            "root_vowel": "/aɪ/ (tense diphthong)",
            "derived_vowel": "/ɪ/ (lax vowel)",
            "explanation": "The tense diphthong in 'derive' becomes a lax vowel in 'derivative'. This is a classic case of trisyllabic laxing.",
            "has_laxed": True
        },
        {
            "word": "serenity",
            "root": "serene",
            "root_vowel": "/iː/ (tense/long vowel)",
            "derived_vowel": "/ɛ/ (lax vowel)",
            "explanation": "The long vowel in 'serene' becomes a lax vowel in 'serenity'. This is a classic case of trisyllabic laxing.",
            "has_laxed": True
        },
        {
            "word": "pleasant",
            "root": "please",
            "root_vowel": "/iː/ (tense/long vowel)",
            "derived_vowel": "/ɛ/ (lax vowel)",
            "explanation": "The long vowel in 'please' becomes a lax vowel in 'pleasant'. Laxing has occurred.",
            "has_laxed": True
        },
        {
            "word": "gratitude",
            "root": "grateful (from Latin grātus)",
            "root_vowel": "/eɪ/ (tense vowel in related 'grateful')",
            "derived_vowel": "/æ/ (lax vowel)",
            "explanation": "Related to 'grateful', which retains a tense vowel, 'gratitude' shows a laxed vowel due to its syllabic structure.",
            "has_laxed": True
        },
        {
            "word": "shadow",
            "root": "Old English 'sceadu'",
            "root_vowel": "/æ/ (lax/short vowel)",
            "derived_vowel": "/æ/ (lax/short vowel)",
            "explanation": "The root vowel in Old English ('sceadu') was already short. Since it did not start with a tense vowel, it could not undergo laxing (shortening).",
            "has_laxed": False
        }
    ]

    # Variable to hold the final answer
    answer = None

    print("Analyzing each word for vowel laxing:")
    print("=" * 45)

    for item in word_analysis:
        print(f"Word: {item['word']}")
        print(f"Root: {item['root']}")
        print(f"Vowel Change: {item['root_vowel']} -> {item['derived_vowel']}")
        print(f"Analysis: {item['explanation']}")
        
        if not item["has_laxed"]:
            answer = item['word']
        
        print("-" * 45)

    # Print the final answer clearly
    if answer:
        print(f"\nConclusion: The word that has not undergone laxing is '{answer}'.")
    else:
        print("\nCould not determine the answer from the analysis.")

find_word_without_laxing()