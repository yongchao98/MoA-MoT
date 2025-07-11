def simulate_primate_response(stimulus_location):
    """
    Simulates the primate's response to a visual stimulus based on a specific brain lesion.

    The lesion is to the right optic radiation, sparing Meyer's loop, which disrupts
    conscious vision in the lower left visual field. Subcortical pathways for
    visually-guided action remain intact.
    """
    print(f"--- Simulating Experiment ---")
    print(f"1. A purple stimulus is presented in the: {stimulus_location.replace('_', ' ')}")

    # Step 2: Determine conscious perception based on the primary visual pathway (LGN to V1)
    lesioned_conscious_pathway = "lower_left_quadrant"
    conscious_report = ""
    print("\n2. Analyzing pathway for conscious vision (to Visual Cortex)...")
    if stimulus_location == lesioned_conscious_pathway:
        print("   - Result: Pathway is ABLATED.")
        conscious_report = "Presses the 'no stimulus present' button."
    else:
        print("   - Result: Pathway is INTACT.")
        conscious_report = "Reports seeing the stimulus."

    # Step 3: Determine motor action based on intact subcortical pathways (e.g., to superior colliculus)
    motor_action = ""
    print("\n3. Analyzing pathway for visually-guided action (Subcortical)...")
    # This pathway is presumed intact for all visual fields.
    print("   - Result: Pathway is INTACT.")
    # A stimulus in the left visual field prompts a reach with the contralateral (left) hand.
    hand_for_reaching = "left" if "left" in stimulus_location else "right"
    motor_action = f"Accurately reaches with its {hand_for_reaching} hand for the target."

    # Step 4: Combine the observations to form a conclusion
    print("\n--- Predicted Behavior ---")
    print(f"Conscious Report: The primate {conscious_report}")
    print(f"Motor Action: The primate {motor_action}")

    print("\n--- Final Conclusion ---")
    if "no stimulus present" in conscious_report and "Accurately reaches" in motor_action:
        conclusion = f"Blindsight for stimuli in the {stimulus_location.replace('_', ' ')} in a non-verbal primate"
        print(f"The combination of accurate reaching without conscious awareness demonstrates: {conclusion}.")
    else:
        conclusion = "A different phenomenon (Normal Vision or Pure Blindness) is demonstrated."
        print(conclusion)


# Run the simulation for the specific scenario in the question.
simulate_primate_response('lower_left_quadrant')
<<<A>>>