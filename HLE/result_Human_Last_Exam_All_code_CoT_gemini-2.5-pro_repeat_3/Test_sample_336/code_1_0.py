def analyze_latin_sentence():
    """
    Analyzes the Latin sentence to find adjectives modifying "mariti".
    """
    sentence = "Sed Maxentium suppositum ferunt arte muliebri tenere mariti animum laborantis auspicio gratissimi partus coepti a puero."
    
    # Grammatical analysis of the target noun
    target_noun = "mariti"
    target_case = "genitive"
    target_number = "singular"
    target_gender = "masculine"
    
    print(f"Analyzing the sentence for adjectives modifying '{target_noun}'.")
    print(f"'{target_noun}' is {target_case} {target_number} {target_gender}, meaning 'of the husband'.\n")

    # List of potential adjectives and participles with their analysis
    adjectives = {
        "muliebri": "Ablative singular. Modifies 'arte' (arte muliebri = by womanly art). Does not modify 'mariti'.",
        "laborantis": "Genitive singular present participle. Agrees with 'mariti'. 'mariti laborantis' means 'of the struggling husband'. This is a match.",
        "gratissimi": "Genitive singular superlative. Modifies 'partus' (gratissimi partus = of a most pleasing birth). Does not modify 'mariti'.",
        "coepti": "Genitive singular perfect participle. Modifies 'partus' (partus coepti = of a birth undertaken). Does not modify 'mariti'."
    }
    
    modifying_adjectives = []
    print("Checking each potential adjective/participle:")
    for adj, analysis in adjectives.items():
        print(f"- {adj}: {analysis}")
        if "This is a match." in analysis:
            modifying_adjectives.append(adj)

    count = len(modifying_adjectives)
    
    print("\n--- Conclusion ---")
    print(f"The number of adjectives modifying '{target_noun}' is: {count}")
    if count > 0:
        print(f"The adjective is: {', '.join(modifying_adjectives)}")

analyze_latin_sentence()