def solve_latin_grammar():
    """
    Analyzes the Latin sentence to find adjectives modifying 'mariti'.
    """
    sentence_structure = {
        'mariti': {'case': 'genitive', 'number': 'singular', 'gender': 'masculine', 'type': 'noun'},
        'muliebri': {'case': 'ablative', 'number': 'singular', 'gender': 'feminine', 'type': 'adjective', 'modifies': 'arte'},
        'laborantis': {'case': 'genitive', 'number': 'singular', 'gender': 'all', 'type': 'participle'},
        'gratissimi': {'case': 'genitive', 'number': 'singular', 'gender': 'masculine/neuter', 'type': 'adjective', 'modifies': 'partus'},
        'coepti': {'case': 'genitive', 'number': 'singular', 'gender': 'masculine/neuter', 'type': 'participle', 'modifies': 'partus'},
    }

    target_noun = "mariti"
    target_properties = sentence_structure[target_noun]

    potential_modifiers = ['muliebri', 'laborantis', 'gratissimi', 'coepti']
    modifying_adjectives = []

    # Check for agreement in case and number, and consider context.
    # 'laborantis' is genitive singular and contextually modifies 'mariti'.
    modifying_adjectives.append('laborantis')
    
    # 'muliebri' is ablative and modifies 'arte'.
    # 'gratissimi' and 'coepti' are genitive singular but modify 'partus' due to proximity and sense.

    count = len(modifying_adjectives)
    
    print(f"The noun 'mariti' is in the genitive singular masculine case.")
    print(f"To modify it, an adjective must also be in the genitive singular masculine case.")
    print(f"Based on grammatical analysis and context:")
    print(f" - 'muliebri' modifies 'arte'.")
    print(f" - 'gratissimi' and 'coepti' modify 'partus'.")
    print(f" - 'laborantis' is the only word that modifies 'mariti'.")
    print("\nResult:")
    print(f"Number of adjectives modifying 'mariti': {count}")
    if modifying_adjectives:
        print(f"List of adjective(s): {', '.join(modifying_adjectives)}")

solve_latin_grammar()