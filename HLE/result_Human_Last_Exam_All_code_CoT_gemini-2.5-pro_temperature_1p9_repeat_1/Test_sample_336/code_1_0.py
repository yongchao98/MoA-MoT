def analyze_mariti_modifiers():
    """
    This script analyzes a Latin sentence to find the number of adjectives
    modifying the noun 'mariti'. It codifies the grammatical analysis step by step.
    """
    sentence = "Sed Maxentium suppositum ferunt arte muliebri tenere mariti animum laborantis auspicio gratissimi partus coepti a puero."
    target_word = "mariti"

    print(f"The task is to find how many adjectives in the sentence modify '{target_word}'.")
    print(f"'{target_word}' is genitive singular masculine ('of the husband'). Modifying adjectives must agree with it.")
    print("\nWe will check each potential adjective and assign 1 if it modifies 'mariti' and 0 otherwise.")

    # Does 'muliebri' modify 'mariti'?
    # No. 'muliebri' modifies 'arte' (ablative case).
    count_muliebri = 0
    print(f"\n- Checking 'muliebri': It is in the ablative phrase 'arte muliebri'. It does not modify the genitive '{target_word}'. Count = {count_muliebri}")

    # Does 'laborantis' modify 'mariti'?
    # Yes. 'laborantis' is a genitive participle that agrees with 'mariti'.
    count_laborantis = 1
    print(f"- Checking 'laborantis': It is a genitive singular participle that modifies '{target_word}' ('of the struggling husband'). Count = {count_laborantis}")

    # Does 'gratissimi' modify 'mariti'?
    # No. 'gratissimi' modifies 'partus' (genitive case).
    count_gratissimi = 0
    print(f"- Checking 'gratissimi': It is a genitive adjective that modifies 'partus' ('of a most pleasing birth'). Count = {count_gratissimi}")

    # Does 'coepti' modify 'mariti'?
    # No. 'coepti' modifies 'partus' (genitive case).
    count_coepti = 0
    print(f"- Checking 'coepti': It is a genitive participle that modifies 'partus' ('of a birth begun'). Count = {count_coepti}")

    # Calculate the total based on the analysis
    total_adjectives = count_muliebri + count_laborantis + count_gratissimi + count_coepti

    # Display the final "equation" as requested
    print("\nTo find the total, we sum the counts for each adjective:")
    print(f"Final Equation: {count_muliebri} + {count_laborantis} + {count_gratissimi} + {count_coepti} = {total_adjectives}")

    print(f"\nConclusion: There is {total_adjectives} adjective modifying 'mariti', which is 'laborantis'.")

analyze_mariti_modifiers()
<<<B>>>