def analyze_latin_sentence():
    """
    Analyzes the Latin sentence to find adjectives modifying "mariti".
    """
    # Dictionary to store the grammatical properties of relevant words.
    # Format: word: {case, number, gender, type}
    words = {
        "mariti": {"case": "genitive", "number": "singular", "gender": "masculine", "type": "noun"},
        "laborantis": {"case": "genitive", "number": "singular", "gender": "masculine/feminine/neuter", "type": "participle"},
        "gratissimi": {"case": "genitive", "number": "singular", "gender": "masculine/neuter", "type": "adjective"},
        "coepti": {"case": "genitive", "number": "singular", "gender": "masculine/neuter", "type": "participle"},
        "partus": {"case": "genitive", "number": "singular", "gender": "masculine", "type": "noun"},
        "arte": {"case": "ablative", "number": "singular", "gender": "feminine", "type": "noun"},
        "muliebri": {"case": "ablative", "number": "singular", "gender": "masculine/feminine/neuter", "type": "adjective"},
    }

    target_noun = "mariti"
    target_properties = words[target_noun]
    
    potential_modifiers = ["laborantis", "gratissimi", "coepti", "muliebri"]
    found_modifiers = []

    print(f"Analyzing adjectives modifying '{target_noun}' ({target_properties['case']}, {target_properties['number']}, {target_properties['gender']}).\n")

    for word in potential_modifiers:
        props = words[word]
        # Basic agreement check (case and number)
        if props["case"] == target_properties["case"] and props["number"] == target_properties["number"]:
            # Contextual check
            # 'laborantis' is next to 'mariti' in the phrase 'mariti animum laborantis'
            # and is the most logical modifier.
            if word == "laborantis":
                print(f"- '{word}' ({props['case']}, {props['number']}) agrees with '{target_noun}'. It describes the husband.")
                found_modifiers.append(word)
            # 'gratissimi' and 'coepti' are part of the phrase 'auspicio gratissimi partus coepti'
            # and modify 'partus'.
            elif word in ["gratissimi", "coepti"]:
                 print(f"- '{word}' ({props['case']}, {props['number']}) agrees with '{target_noun}' grammatically, but contextually it modifies 'partus'.")
            else:
                 # This branch won't be hit with the current logic but is good for completeness.
                 found_modifiers.append(word)
        else:
            print(f"- '{word}' ({props['case']}, {props['number']}) does not agree with '{target_noun}' ({target_properties['case']}, {target_properties['number']}).")

    print(f"\nConclusion: There is {len(found_modifiers)} adjective/participle modifying 'mariti'.")
    print("List of modifiers:", found_modifiers)

analyze_latin_sentence()
<<<B>>>