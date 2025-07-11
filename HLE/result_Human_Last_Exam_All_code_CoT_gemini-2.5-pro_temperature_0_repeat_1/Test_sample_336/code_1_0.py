def solve_latin_grammar():
    """
    Analyzes the Latin sentence to find adjectives modifying 'mariti'.
    """
    sentence = "Sed Maxentium suppositum ferunt arte muliebri tenere mariti animum laborantis auspicio gratissimi partus coepti a puero."
    target_noun = "mariti"
    target_noun_grammar = {
        "case": "Genitive",
        "number": "Singular",
        "gender": "Masculine"
    }

    adjectives_analysis = {
        "muliebri": "Ablative Singular, modifies 'arte' (art). Does not modify 'mariti'.",
        "laborantis": "Genitive Singular, agrees with 'mariti'. Modifies 'mariti' (husband).",
        "gratissimi": "Genitive Singular, agrees with 'partus' (birth). Modifies 'partus'.",
        "coepti": "Genitive Singular, agrees with 'partus' (birth). Modifies 'partus'."
    }

    modifying_adjectives = []
    for adj, analysis in adjectives_analysis.items():
        if "Modifies 'mariti'" in analysis:
            modifying_adjectives.append(adj)

    count = len(modifying_adjectives)
    
    print(f"The noun in question is 'mariti', which is genitive, singular, and masculine.")
    print(f"We need to find adjectives that agree with it in case, number, and gender.")
    print("\nAnalysis of potential adjectives:")
    print(f"- 'muliebri' is ablative singular and modifies 'arte'.")
    print(f"- 'laborantis' is genitive singular and modifies 'mariti'.")
    print(f"- 'gratissimi' is genitive singular and modifies 'partus'.")
    print(f"- 'coepti' is genitive singular and modifies 'partus'.")
    
    print(f"\nBased on the analysis, there is only {count} adjective that modifies 'mariti'.")
    print(f"The adjective is: {', '.join(modifying_adjectives)}")

solve_latin_grammar()