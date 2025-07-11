def analyze_latin_sentence():
    """
    Analyzes the Latin sentence to find adjectives modifying 'mariti'.
    """
    # The noun to be analyzed and its grammatical properties
    target_noun = {
        "word": "mariti",
        "case": "genitive",
        "number": "singular",
        "gender": "masculine"
    }

    # List of potential adjectives/participles and their analysis
    potential_modifiers = [
        {"word": "muliebri", "modifies": "arte", "reason": "Agrees with 'arte' (ablative, singular, feminine)."},
        {"word": "laborantis", "modifies": "mariti", "reason": "Agrees with 'mariti' (genitive, singular, masculine)."},
        {"word": "gratissimi", "modifies": "partus", "reason": "Agrees with 'partus' (genitive, singular, masculine)."},
        {"word": "coepti", "modifies": "partus", "reason": "Agrees with 'partus' (genitive, singular, masculine)."}
    ]

    modifying_adjectives = []
    print("Analyzing which adjectives modify 'mariti' (genitive, singular, masculine):")
    print("-" * 60)

    for item in potential_modifiers:
        print(f"Checking '{item['word']}':")
        if item['modifies'] == target_noun['word']:
            print(f"  -> YES, it modifies '{target_noun['word']}'.")
            print(f"  -> Reason: {item['reason']}")
            modifying_adjectives.append(item['word'])
        else:
            print(f"  -> NO, it modifies '{item['modifies']}'.")
            print(f"  -> Reason: {item['reason']}")
        print()

    count = len(modifying_adjectives)
    print("-" * 60)
    print(f"Final Count: There is {count} adjective/participle modifying '{target_noun['word']}'.")
    if count > 0:
        print(f"The modifying word is: {', '.join(modifying_adjectives)}")

analyze_latin_sentence()
<<<B>>>