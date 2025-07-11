def simulate_primate_response(lesion_location, stimulus_quadrant):
    """
    Simulates a primate's response to a visual stimulus based on a brain lesion.

    Args:
        lesion_location (str): A description of the brain lesion.
        stimulus_quadrant (str): The quadrant where the stimulus is presented.
    """
    # Step 1: Determine the affected visual field based on neuroanatomy.
    # A lesion on the right optic radiation affects the left visual field.
    # Sparing Meyer's loop means the fibers for the INFERIOR field are damaged.
    affected_quadrant = "lower left quadrant"

    print(f"Condition: The primate has a lesion in the '{lesion_location}'.")
    print(f"Resulting Deficit: This damages the primary visual pathway for the {affected_quadrant}.")
    print("-" * 30)
    print(f"Experiment: A purple stimulus is presented in the {stimulus_quadrant}.")
    print("-" * 30)

    # Step 2: Check if the stimulus is in the cortically "blind" field.
    if stimulus_quadrant == affected_quadrant:
        # Step 3: Describe the two conflicting behaviors seen in blindsight.
        conscious_report = "Presses the 'no trial' tone, indicating no conscious perception."
        motor_action = "Accurately reaches with its left hand for the stimulus."

        print("Observed Behaviors:")
        print(f"1. Conscious Report: {conscious_report}")
        print(f"2. Motor Action: {motor_action}")
        print("-" * 30)

        # Step 4: State the final conclusion.
        conclusion = "Blindsight for stimuli in the lower left quadrant in a non-verbal primate"
        print(f"Conclusion: The primate demonstrates: {conclusion}.")
    else:
        print("Conclusion: The primate would demonstrate normal vision in this quadrant.")


# --- Main Execution ---
# Define the scenario from the problem description.
primate_lesion = "right optic radiation, sparing Meyer's loop"
stimulus_presentation = "lower left quadrant"

simulate_primate_response(primate_lesion, stimulus_presentation)