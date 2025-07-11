def analyze_latin_sentence():
    """
    Analyzes a Latin sentence to find which adjectives modify the noun "mariti".
    """
    sentence = "Sed Maxentium suppositum ferunt arte muliebri tenere mariti animum laborantis auspicio gratissimi partus coepti a puero."
    target_noun = "mariti"
    
    analysis = []
    
    # Step 1: Analyze the target noun 'mariti'
    analysis.append(f"Analyzing the sentence: \"{sentence}\"")
    analysis.append(f"The noun in question is '{target_noun}'.")
    analysis.append(f"Grammatical analysis of '{target_noun}':")
    analysis.append("- Case: Genitive (of the husband)")
    analysis.append("- Number: Singular")
    analysis.append("- Gender: Masculine")
    analysis.append("An adjective must agree with 'mariti' in case, number, and gender.")
    analysis.append("\n--- Checking potential adjectives/participles ---")

    adjectives = {
        "muliebri": "Ablative, singular, feminine. Modifies 'arte' (arte muliebri - 'by a womanly artifice'). Does NOT modify 'mariti'.",
        "laborantis": "Present participle, genitive, singular, all genders. It AGREES with 'mariti'. The phrase 'animum mariti laborantis' means 'the mind/affections of the struggling husband'. The meaning and grammar align perfectly.",
        "gratissimi": "Superlative adjective, genitive, singular, masculine/neuter. It could technically agree with 'mariti'. However, it is positioned with 'partus' (of the birth). The phrase is 'auspicio gratissimi partus' ('by the good omen of a most pleasing birth'). It clearly modifies 'partus'. Does NOT modify 'mariti'.",
        "coepti": "Perfect passive participle, genitive, singular, masculine/neuter. It could also technically agree with 'mariti'. However, it follows 'partus' and forms the phrase 'partus coepti a puero' ('a birth begun by a boy'). It modifies 'partus'. Does NOT modify 'mariti'."
    }
    
    modifying_adjectives = []
    for adj, desc in adjectives.items():
        analysis.append(f"\n- Adjective/Participle: '{adj}'")
        analysis.append(f"  Analysis: {desc}")
        if "AGREES with 'mariti'" in desc or "modifies 'mariti'" in desc:
            modifying_adjectives.append(adj)

    analysis.append("\n--- Conclusion ---")
    analysis.append(f"The total number of adjectives modifying 'mariti' is: {len(modifying_adjectives)}")
    
    if modifying_adjectives:
        analysis.append("The adjective modifying 'mariti' is:")
        for adj in modifying_adjectives:
            analysis.append(f"- {adj}")
    else:
        analysis.append("There are no adjectives modifying 'mariti'.")

    # Print the analysis step-by-step
    for line in analysis:
        print(line)

analyze_latin_sentence()
<<<B>>>