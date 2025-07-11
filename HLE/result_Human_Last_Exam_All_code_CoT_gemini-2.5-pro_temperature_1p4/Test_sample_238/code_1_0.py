def explain_guarani_effected_objects():
    """
    Explains and demonstrates the interaction between Guarani's nominal tense
    and effected objects in sentences.
    """
    # 1. Define the linguistic concepts and markers
    effected_object_def = "An object whose existence is a result of the verb's action."
    destinative_marker = "-rã"
    destinative_meaning = "Marks a noun as 'future', 'potential', or 'destined-to-be'."
    
    # 2. Explain the grammatical rule
    print("Step 1: Understanding the Rule")
    print(f"In Guarani, an '{effected_object_def}' must be grammatically marked to show its status.")
    print("If the action has not yet happened, the object does not yet exist. It is a potential object.")
    print(f"Therefore, it must be marked with the destinative suffix '{destinative_marker}', which means '{destinative_meaning}'.\n")
    
    # 3. Provide a concrete example
    print("Step 2: Building an Example")
    base_noun = "óga"
    meaning_noun = "house"
    verb_phrase = "a-japo-ta"
    meaning_verb = "I will make/build"
    
    print(f"Let's say we want to express 'I will build a house'.")
    print(f"The base noun for 'house' is '{base_noun}'.")
    print("Since the house will be 'effected' (created) by the future action, we must modify the noun.\n")
    
    # 4. Construct the final "equation" showing each part
    print("Step 3: The Final Equation")
    print("To form the correct noun phrase for the 'house-to-be-built', we apply the rule:")
    
    # Print each component of the final grammatical construction
    final_form = base_noun + destinative_marker
    print(f"Base Noun ('{base_noun}') + Destinative Suffix ('{destinative_marker}') => Final Form ('{final_form}')")
    
    print(f"\nThe full sentence would be: '{verb_phrase} peteĩ {final_form}'.")

# Execute the explanation
explain_guarani_effected_objects()