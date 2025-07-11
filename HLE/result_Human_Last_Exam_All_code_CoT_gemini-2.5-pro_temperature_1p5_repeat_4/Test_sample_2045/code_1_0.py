def diagnose_visual_field_defect(lesion_hemisphere, lesion_pathway, behavior_description):
    """
    Diagnoses a visual field defect based on lesion location and describes the resulting phenomenon.

    Args:
        lesion_hemisphere (str): The hemisphere of the lesion ('right' or 'left').
        lesion_pathway (str): The part of the optic radiation affected ('meyers_loop' or 'parietal_pathway').
        behavior_description (dict): A dictionary describing the subject's behavior.
    """

    print("Step 1: Analyzing the lesion location...")
    # Determine contralateral visual field
    if lesion_hemisphere == 'right':
        affected_field_side = 'left'
    else:
        affected_field_side = 'right'
    print(f"A lesion in the {lesion_hemisphere} hemisphere affects the contralateral (opposite) visual field: the '{affected_field_side}' visual field.")

    # Determine vertical quadrant
    if lesion_pathway == 'parietal_pathway':
        # This pathway (outside Meyer's loop) serves the lower visual field.
        affected_field_vertical = 'lower'
    else: # meyers_loop
        # Meyer's loop serves the upper visual field.
        affected_field_vertical = 'upper'
    print(f"The lesion is in the {lesion_pathway}, which corresponds to the '{affected_field_vertical}' visual field.")

    print("\nStep 2: Combining the findings to identify the blind quadrant...")
    affected_quadrant = f"{affected_field_vertical} {affected_field_side}"
    print(f"Therefore, the lesion results in a visual defect in the '{affected_quadrant}' quadrant.")

    print("\nStep 3: Interpreting the observed behavior...")
    conscious_report = behavior_description['conscious_report']
    motor_action = behavior_description['motor_action']
    print(f"The primate reports '{conscious_report}' when stimuli are in the affected quadrant.")
    print(f"However, the primate is able to '{motor_action}' in the same quadrant.")
    print("This disconnect between conscious perception and accurate physical interaction is known as Blindsight.")

    print("\n--- FINAL CONCLUSION ---")
    print(f"The primate demonstrates: Blindsight for stimuli in the {affected_quadrant} quadrant in a non-verbal primate.")


# Parameters from the user's problem description
lesion_in_hemisphere = 'right'
# The lesion is "outside the Meyer's loop", which is the parietal pathway
lesion_in_pathway = 'parietal_pathway'
primate_behavior = {
    "conscious_report": "no stimulus present (by pressing a specific button)",
    "motor_action": "accurately reach for the target with its left hand"
}

diagnose_visual_field_defect(lesion_in_hemisphere, lesion_in_pathway, primate_behavior)
