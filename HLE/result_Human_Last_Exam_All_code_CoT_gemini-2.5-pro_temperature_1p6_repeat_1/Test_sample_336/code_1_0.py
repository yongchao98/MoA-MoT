def analyze_latin_sentence():
    """
    Analyzes a Latin sentence to find adjectives modifying a specific noun.
    """
    # The noun we are interested in is "mariti".
    # Grammatical analysis shows "mariti" is genitive, singular, masculine.
    target_noun = {
        'word': 'mariti',
        'case': 'genitive',
        'number': 'singular',
        'gender': 'masculine'
    }

    # Let's check the potential adjectives in the sentence.
    # An adjective must agree with the noun it modifies in case, number, and gender.
    adjectives_to_check = [
        {
            'word': 'laborantis',
            # This form is genitive singular for all genders.
            'case': 'genitive',
            'number': 'singular',
            'gender': ['masculine', 'feminine', 'neuter'],
            'meaning': 'struggling'
        },
        {
            'word': 'gratissimi',
            # This form is genitive singular for masculine/neuter.
            'case': 'genitive',
            'number': 'singular',
            'gender': ['masculine', 'neuter'],
            'meaning': 'most welcome'
        },
        {
            'word': 'coepti',
            # This form is genitive singular for masculine/neuter.
            'case': 'genitive',
            'number': 'singular',
            'gender': ['masculine', 'neuter'],
            'meaning': 'begun/undertaken'
        }
    ]

    # Another genitive noun in the sentence is "partus" (of a birth).
    # It is also genitive, singular, masculine.
    other_noun = {
        'word': 'partus',
        'case': 'genitive',
        'number': 'singular',
        'gender': 'masculine'
    }

    modifying_adjectives = []

    print("Analyzing which adjectives modify 'mariti' (genitive, singular, masculine):")
    print("-" * 60)

    # Check each adjective
    for adj in adjectives_to_check:
        # Check for grammatical agreement with 'mariti'
        if (adj['case'] == target_noun['case'] and
            adj['number'] == target_noun['number'] and
            target_noun['gender'] in adj['gender']):
            
            # Now, consider context and meaning (semantic check)
            if adj['word'] == 'laborantis':
                print(f"'{adj['word']}': Grammatically agrees. Semantically, 'mariti laborantis' ('of the struggling husband') makes sense.")
                modifying_adjectives.append(adj['word'])
            elif adj['word'] == 'gratissimi':
                print(f"'{adj['word']}': Grammatically agrees. Semantically, it is positioned next to 'partus'. 'gratissimi partus' ('of a most welcome birth') is a much better fit.")
            elif adj['word'] == 'coepti':
                print(f"'{adj['word']}': Grammatically agrees. Semantically, it modifies 'partus'. 'partus coepti' ('birth undertaken') makes sense, while 'mariti coepti' ('husband undertaken') does not.")
        else:
            # This case won't be hit with our simplified data, but is good practice.
            print(f"'{adj['word']}' does not grammatically agree with 'mariti'.")

    print("-" * 60)
    print(f"Total number of adjectives modifying 'mariti': {len(modifying_adjectives)}")
    if modifying_adjectives:
        print(f"The adjective is: {', '.join(modifying_adjectives)}")
    else:
        print("There are no adjectives modifying 'mariti'.")

analyze_latin_sentence()
<<<B>>>