import sys

def analyze_primate_case():
    """
    Analyzes the neurological case of the primate to determine the demonstrated phenomenon.
    """
    # 1. Define the anatomical and behavioral parameters from the prompt.
    lesion_hemisphere = "Right"
    lesion_pathway = "Optic Radiation"
    spared_structure = "Meyer's Loop"
    
    stimulus_quadrant_x = "Left"
    stimulus_quadrant_y = "Lower"
    
    # Behavior 1: Motor response to the stimulus
    motor_response_action = "Accurate reach with left hand"
    motor_response_indicates = "Preserved subconscious visual processing for motor guidance"
    
    # Behavior 2: Conscious report about the stimulus
    conscious_report_action = "Presses 'no trial' button"
    conscious_report_indicates = "Lack of conscious visual perception"

    # 2. Deduce the affected visual field from the lesion location.
    # Visual pathways are contralateral (opposite side).
    if lesion_hemisphere == "Right":
        affected_visual_field_x = "Left"
    else:
        affected_visual_field_x = "Right"
        
    # Meyer's loop carries SUPERIOR visual field data. The lesion is OUTSIDE it,
    # affecting the fibers that carry INFERIOR visual field data.
    if spared_structure == "Meyer's Loop":
        affected_visual_field_y = "Lower"
    else:
        # This is a hypothetical for a different lesion
        affected_visual_field_y = "Upper"

    affected_quadrant = f"{affected_visual_field_y} {affected_visual_field_x} Quadrant"

    # 3. Print the step-by-step logical deduction.
    print("--- Case Analysis ---")
    print(f"Step 1: The lesion is in the {lesion_hemisphere} Hemisphere's {lesion_pathway}, sparing {spared_structure}.")
    print(f"Step 2: A right-sided brain lesion affects the left visual field.")
    print(f"Step 3: Sparing Meyer's Loop means the pathways for the lower visual field are damaged.")
    print(f"Step 4: Therefore, the predicted area of conscious blindness is the '{affected_quadrant}'.")
    print("-" * 25)
    print("--- Behavioral Observation ---")
    print(f"Stimulus Location: A stimulus is presented in the '{stimulus_quadrant_y} {stimulus_quadrant_x} Quadrant'.")
    print(f"Observation 1 (Motor): The primate shows an '{motor_response_action}'. This suggests: {motor_response_indicates}.")
    print(f"Observation 2 (Report): The primate also performs the action: '{conscious_report_action}'. This suggests: {conscious_report_indicates}.")
    print("-" * 25)
    print("--- Conclusion ---")
    print("The combination of accurate action without conscious awareness in a blind field is known as Blindsight.")
    print(f"Final Diagnosis: The primate demonstrates Blindsight for stimuli in the {affected_quadrant}.")

if __name__ == '__main__':
    analyze_primate_case()