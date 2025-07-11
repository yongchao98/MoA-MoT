def analyze_latin_adjectives():
    """
    Analyzes a Latin sentence to find which adjectives modify the noun "mariti".
    """
    target_noun = {
        "word": "mariti",
        "case": "genitive",
        "number": "singular",
        "gender": "masculine"
    }

    adjectives_to_check = {
        "muliebri": {"case": "ablative", "number": "singular", "modifies": "arte"},
        "laborantis": {"case": "genitive", "number": "singular", "modifies": "mariti"},
        "gratissimi": {"case": "genitive", "number": "singular", "modifies": "partus"},
        "coepti": {"case": "genitive", "number": "singular", "modifies": "partus"}
    }

    modifying_list = []

    print("--- Grammatical Analysis ---")
    print(f"The goal is to find adjectives modifying '{target_noun['word']}'.")
    print(f"'{target_noun['word']}' is {target_noun['case']}, {target_noun['number']}, and {target_noun['gender']}.")
    print("Any modifying adjective must agree in case, number, and gender.\n")

    for adj, properties in adjectives_to_check.items():
        if (properties["case"] == target_noun["case"] and
            properties["number"] == target_noun["number"] and
            properties["modifies"] == target_noun["word"]):
            print(f"- '{adj}' is {properties['case']} {properties['number']}. It agrees with and modifies '{target_noun['word']}'.")
            modifying_list.append(adj)
        elif properties["modifies"] != target_noun["word"]:
            print(f"- '{adj}' is {properties['case']} {properties['number']}. It does not modify '{target_noun['word']}' because it modifies '{properties['modifies']}'.")

    count = len(modifying_list)

    print("\n--- Final Answer ---")
    print(f"The number of adjectives modifying 'mariti' is: {count}")
    if count > 0:
        print(f"The adjective is: {', '.join(modifying_list)}")

analyze_latin_adjectives()
<<<B>>>