def sound_hierarchy_for_dance_scene():
    """
    This function defines and prints the hierarchy of sound importance in a typical film dance scene.
    """
    # The four primary elements of film sound. "Ambience" is added to complete the user's list.
    # The list is ordered from most to least important for a dance scene.
    ranked_elements = ["Music", "Sound effects", "Speech", "Ambience"]

    # The prompt asks to "output each number in the final equation".
    # We will format the output as a visual equation showing the ranking.
    equation_parts = []
    for i, element in enumerate(ranked_elements, 1):
        # We pair each element with its rank number.
        equation_parts.append(f"{i}. {element}")

    # Join the parts with a '+' to create the visual "equation".
    final_equation = " + ".join(equation_parts)

    print("In a film's dance scene, the hierarchy of sound importance can be represented by the following equation:")
    print(final_equation)

sound_hierarchy_for_dance_scene()