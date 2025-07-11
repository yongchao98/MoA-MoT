def analyze_latin_sentence():
    """
    Analyzes a Latin sentence to find and count the adjectives modifying "mariti".
    """
    # The noun "mariti" is genitive, singular, masculine.
    # We will check which potential adjectives agree with it.

    candidates = {
        "muliebri": {
            "modifies": False,
            "reason": "Does not modify 'mariti'. It is ablative singular and modifies 'arte'."
        },
        "laborantis": {
            "modifies": True,
            "reason": "Modifies 'mariti'. It is a genitive singular participle, agreeing with 'mariti'."
        },
        "gratissimi": {
            "modifies": False,
            "reason": "Does not modify 'mariti'. It is genitive singular and modifies 'partus'."
        },
        "coepti": {
            "modifies": False,
            "reason": "Does not modify 'mariti'. It is genitive singular and modifies 'partus'."
        }
    }

    print("To find the number of adjectives modifying 'mariti', we check each potential candidate:")
    
    count = 0
    equation_parts = []
    
    for word, analysis in candidates.items():
        if analysis["modifies"]:
            value = 1
        else:
            value = 0
        
        print(f"- Checking '{word}': {analysis['reason']} This adds {value} to the total.")
        count += value
        equation_parts.append(str(value))

    # The user requested to see the numbers in the final equation.
    # The order corresponds to the analysis order: muliebri, laborantis, gratissimi, coepti.
    final_equation = " + ".join(equation_parts)

    print("\nFinal Count Calculation:")
    print(f"{final_equation} = {count}")
    
    if count == 1:
        modifying_adjective = [word for word, analysis in candidates.items() if analysis["modifies"]][0]
        print(f"\nThere is 1 adjective ('{modifying_adjective}') that modifies 'mariti'.")
    else:
        print(f"\nThere are {count} adjectives that modify 'mariti'.")

analyze_latin_sentence()