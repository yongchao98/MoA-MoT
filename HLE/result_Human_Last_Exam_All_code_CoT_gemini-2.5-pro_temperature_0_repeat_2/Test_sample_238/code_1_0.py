def apply_guarani_nominal_marker(noun, object_type):
    """
    This function demonstrates the grammatical rule for applying nominal tense
    markers to nouns in Guarani based on whether they are effected objects.

    An 'effected object' is created by the verb's action (e.g., building a house).
    """
    marker = ""
    explanation = ""

    if object_type == "effected":
        # Effected objects are the future result of an action.
        # They take the destinative/future marker '-rã'.
        marker = "-rã"
        explanation = f"The noun '{noun}' is an effected object, so it takes the destinative marker '{marker}'."
    else:
        # Other object types would follow different rules.
        explanation = f"The noun '{noun}' is not an effected object in this context."

    final_form = noun + marker

    print("--- Guarani Grammar Demonstration ---")
    print(f"Base Noun: {noun}")
    print(f"Object Type: {object_type}")
    print(f"Rule Applied: {explanation}")
    print(f"Final Form: {final_form}")


# --- Example ---
# Let's model the object in the phrase "to build a house" (óga apo).
# The "house" (óga) is created by the action, so it's an "effected object".
noun_to_mark = "óga"
type_of_object = "effected"

apply_guarani_nominal_marker(noun_to_mark, type_of_object)