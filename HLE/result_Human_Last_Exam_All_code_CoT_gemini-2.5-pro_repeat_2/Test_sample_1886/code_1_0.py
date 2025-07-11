import textwrap

def analyze_circumflex_functions():
    """
    Analyzes the attested functions of the circumflex accent in French
    to determine which of the given options is not a valid function.
    """
    options = {
        'A': {
            "text": "To indicate contrastive distinctions between closed and open vowel sounds.",
            "is_function": True,
            "explanation": "This is a real phonetic function. For example, 'ê' (as in 'fête') is an open vowel [ɛ], reflecting a historical long vowel."
        },
        'B': {
            "text": "To distinguish words that are pronounced the same way, but have more than one meaning.",
            "is_function": True,
            "explanation": "This is a key grammatical function to distinguish homophones, e.g., 'sur' (on) vs. 'sûr' (sure), or 'du' (of the) vs. 'dû' (due)."
        },
        'C': {
            "text": "To indicate a vowel pronounced as [o] within words originating in Classical Latin.",
            "is_function": True,
            "explanation": "This is an indirect but attested occurrence. The circumflex on 'ô' can mark a historical long vowel (e.g., 'cône' from Latin 'conus') or a lost 's' (e.g., 'côte' from 'costa'). In both cases, it indicates an [o] sound in a word of Latin origin."
        },
        'D': {
            "text": "To distinguish between short and long vowel sounds.",
            "is_function": True,
            "explanation": "This was a primary historical function of the circumflex, though the length distinction is lost in many modern French dialects. Example: 'pâte' (long a) vs. 'patte' (short a)."
        },
        'E': {
            "text": "None of the options.",
            "is_function": None,
            "explanation": "This would be true if all other options were valid functions."
        },
        'F': {
            "text": "To make a word appear more prestigious.",
            "is_function": False,
            "explanation": "This is a sociolinguistic effect, not a prescribed orthographic function. Diacritics are added for linguistic reasons (phonetic, etymological, grammatical), not to confer prestige. Prestige may be a byproduct of 'correct' usage, but it is not the purpose of the rule itself."
        },
        'G': {
            "text": "To indicate where a diphthong has been reduced to a single vowel sound.",
            "is_function": True,
            "explanation": "This is a valid etymological path. The reduction of a diphthong often resulted in a long vowel, which was then marked with a circumflex (e.g., the evolution of the vowel in 'maître' from Old French 'maistre')."
        },
        'H': {
            "text": "To indicate where a sibilant once existed in both the spoken and written language.",
            "is_function": True,
            "explanation": "This is one of the most common and well-known functions, marking a lost 's'. Examples: 'forêt' (from 'forest'), 'château' (from 'chastel')."
        },
        'I': {
            "text": "To indicate where two consecutive vowels in hiatus have been reduced to a single vowel sound.",
            "is_function": True,
            "explanation": "This is a valid etymological function. The classic example is 'âge' which evolved from the Old French 'aage'."
        }
    }

    correct_answer_key = None
    
    print("Analysis of Proposed Functions of the French Circumflex Accent")
    print("="*60)

    # Iterating through sorted keys for consistent order
    for key in sorted(options.keys()):
        option_data = options[key]
        if option_data["is_function"] is None:
            continue

        status = "VALID" if option_data["is_function"] else "INVALID"
        print(f"Option {key}: {status}")
        print(f"  Statement: \"{option_data['text']}\"")
        
        # Wrap the explanation text for better readability
        wrapped_explanation = textwrap.fill(f"  Analysis: {option_data['explanation']}", width=70, subsequent_indent='  ')
        print(wrapped_explanation)
        print("-" * 60)
        
        if not option_data["is_function"]:
            correct_answer_key = key

    print("\nConclusion:")
    if correct_answer_key:
        print(f"The option that has never been an attested orthographic function is '{correct_answer_key}'.")
        print("The purpose of orthographic rules is linguistic (representing sound, history, or grammar), not sociolinguistic (conferring prestige).")
    else:
        print("Error: Could not determine the invalid function from the provided options.")


analyze_circumflex_functions()