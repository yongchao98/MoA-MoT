def diagnose_condition():
    """
    Analyzes neuroanatomical and behavioral data to diagnose a primate's condition.
    """
    # Step 1: Define anatomical and behavioral parameters
    lesion_hemisphere = "right"
    # The part of the optic radiation outside Meyer's loop carries inferior visual field info
    lesion_pathway_component = "parietal (non-Meyer's loop)"
    primate_action = "accurate reaching"
    primate_report = "no stimulus perceived"

    # Assign numerical codes to represent concepts for the "equation"
    RIGHT_HEMISPHERE_CODE = 1
    INFERIOR_FIELD_PATHWAY_CODE = 2
    ACCURATE_ACTION_CODE = 10
    NO_AWARENESS_CODE = 20

    # Step 2: Deduce the affected visual field from the lesion location
    if lesion_hemisphere == "right":
        affected_field_horizontal = "left"
    else:
        affected_field_horizontal = "right"

    if lesion_pathway_component == "parietal (non-Meyer's loop)":
        affected_field_vertical = "lower"
    else:
        affected_field_vertical = "upper"

    affected_quadrant = f"{affected_field_vertical} {affected_field_horizontal} quadrant"

    # Step 3: Deduce the condition from the behavior
    # Blindsight is defined by accurate action without conscious awareness
    if primate_action == "accurate reaching" and primate_report == "no stimulus perceived":
        diagnosis = "Blindsight"
    else:
        diagnosis = "Another condition (e.g., pure blindness or full sight)"

    # Step 4: Print the analysis and the "equation"
    print("--- Analysis Breakdown ---")
    print(f"1. Lesion in the '{lesion_hemisphere}' hemisphere impacts the '{affected_field_horizontal}' visual field.")
    print(f"2. Damage to the '{lesion_pathway_component}' component of the optic radiation impacts the '{affected_field_vertical}' visual field.")
    print(f"--> Resulting deficit is in the: {affected_quadrant}\n")

    print(f"3. Primate demonstrates '{primate_action}' (Action is present).")
    print(f"4. Primate reports '{primate_report}' (Conscious awareness is absent).")
    print("--> This dissociation between action and awareness is characteristic of blindsight.\n")

    print("--- Conceptual Equation ---")
    print("This demonstrates how different factors combine for the final diagnosis:")
    print(f"Location Factors: Hemisphere Code ({RIGHT_HEMISPHERE_CODE}) + Pathway Code ({INFERIOR_FIELD_PATHWAY_CODE}) -> Deficit in Lower Left Quadrant")
    print(f"Behavioral Factors: Action Code ({ACCURATE_ACTION_CODE}) + Awareness Code ({NO_AWARENESS_CODE}) -> Blindsight Behavior")
    print("\n--- Final Conclusion ---")
    print(f"The primate will demonstrate: {diagnosis} for stimuli in the {affected_quadrant}.")

# Run the diagnosis
diagnose_condition()