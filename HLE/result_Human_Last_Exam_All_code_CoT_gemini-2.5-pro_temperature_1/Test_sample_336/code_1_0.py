def analyze_latin_sentence():
    """
    Analyzes a Latin sentence to find which adjectives modify the noun 'mariti'.
    """
    sentence = "Sed Maxentium suppositum ferunt arte muliebri tenere mariti animum laborantis auspicio gratissimi partus coepti a puero."
    target_noun = "mariti"
    candidates = ["muliebri", "laborantis", "gratissimi", "coepti"]

    analysis_results = {
        "mariti": "genitive singular masculine, meaning 'of the husband'. It modifies 'animum' (the mind).",
        "muliebri": "ablative singular feminine. It agrees with and modifies 'arte' (trick). It cannot modify the genitive 'mariti'.",
        "laborantis": "genitive singular participle. It agrees with 'mariti' in case and number. The phrase 'mariti...laborantis' means 'of the struggling husband'. This is a correct match.",
        "gratissimi": "genitive singular superlative adjective. In the phrase 'gratissimi partus' (of a most pleasing birth), it agrees with and modifies 'partus'. It does not modify 'mariti'.",
        "coepti": "genitive singular participle. In the phrase 'partus coepti a puero' (of the birth begun by the boy), it agrees with and modifies 'partus'. It does not modify 'mariti'."
    }

    modifying_adjectives = []

    # Based on grammatical analysis, find the adjective(s) modifying 'mariti'
    if "laborantis" in candidates:
        # Checking the logic described in the analysis dictionary
        # 'laborantis' is genitive singular and logically modifies 'mariti'
        modifying_adjectives.append("laborantis")

    print(f"Analyzing adjectives modifying '{target_noun}':\n")
    print(f"1. The noun '{target_noun}' is {analysis_results[target_noun]}")
    print(f"2. Candidate 'muliebri' is {analysis_results['muliebri']}")
    print(f"3. Candidate 'laborantis' is {analysis_results['laborantis']}")
    print(f"4. Candidate 'gratissimi' is {analysis_results['gratissimi']}")
    print(f"5. Candidate 'coepti' is {analysis_results['coepti']}")
    print("\n--------------------")
    print("Conclusion:")
    
    count = len(modifying_adjectives)
    
    print(f"There is {count} adjective modifying 'mariti'.")
    if modifying_adjectives:
        print(f"The adjective is: {', '.join(modifying_adjectives)}")

analyze_latin_sentence()
<<<B>>>