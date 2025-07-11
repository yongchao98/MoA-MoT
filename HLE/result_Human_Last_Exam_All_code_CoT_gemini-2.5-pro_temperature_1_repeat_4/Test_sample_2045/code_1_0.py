def simulate_primate_response(stimulus_location, stimulus_color):
    """
    Simulates the behavioral response of a primate with a lesion causing
    blindsight in the lower left visual quadrant.

    Args:
        stimulus_location (str): The quadrant where the stimulus appears.
        stimulus_color (str): The color of the stimulus.
    """

    # The lesion to the right superior optic radiation causes a deficit
    # in the lower left visual field.
    affected_quadrant = "lower left"

    print(f"A {stimulus_color} stimulus is given in the {stimulus_location} quadrant.")
    print("---")

    if stimulus_location == affected_quadrant:
        # This is the blindsight condition
        motor_response = "The primate accurately reaches with its left hand for the target."
        awareness_response = "The primate presses the button to indicate 'no stimulus' is present."
        conclusion = "This demonstrates blindsight for stimuli in the lower left quadrant."

    elif stimulus_location in ["upper left", "lower right", "upper right"]:
        # This represents normal vision in unaffected quadrants
        motor_response = "The primate accurately reaches for the target."
        awareness_response = "The primate does NOT press the button, correctly indicating a stimulus is present."
        conclusion = "This demonstrates normal vision."

    else:
        # This is the baseline condition with no stimulus
        motor_response = "The primate does not reach."
        awareness_response = "The primate correctly presses the button to indicate 'no stimulus' is present."
        conclusion = "This demonstrates correct baseline behavior."

    print(f"Predicted Motor Response: {motor_response}")
    print(f"Predicted Awareness Report: {awareness_response}")
    print(f"Conclusion: {conclusion}")


# Run the simulation for the specific scenario described in the problem.
stimulus_location_in_scenario = 'lower left'
stimulus_color_in_scenario = 'purple'
simulate_primate_response(stimulus_location_in_scenario, stimulus_color_in_scenario)