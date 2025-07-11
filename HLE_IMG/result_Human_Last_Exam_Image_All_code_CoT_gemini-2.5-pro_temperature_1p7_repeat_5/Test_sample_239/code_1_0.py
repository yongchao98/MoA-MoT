def quantum_eraser_logic():
    """
    Analyzes the delayed-choice quantum eraser experiment to determine
    the outcome at detector D0.
    """

    print("Analyzing the paths for each detector to determine if 'which-path' information is available or erased.\n")

    # Detector D1 Analysis
    print("Detector D1:")
    print(" - A photon detected at D1 must have followed the 'red' path.")
    print(" - This provides 'which-path' information.")
    print(" - Conclusion for D0 (correlated with D1): No interference pattern.\n")

    # Detector D2 Analysis
    print("Detector D2:")
    print(" - A photon detected at D2 must have followed the 'cyan' path.")
    print(" - This also provides 'which-path' information.")
    print(" - Conclusion for D0 (correlated with D2): No interference pattern.\n")

    # Detectors D3 and D4 Analysis
    print("Detectors D3 and D4:")
    print(" - The 'red' and 'cyan' paths are recombined at beam-splitter BSc before reaching D3 and D4.")
    print(" - A photon detected at D3 or D4 could have come from either path.")
    print(" - 'Which-path' information is erased.")
    print(" - Conclusion for D0 (correlated with D3 or D4): An interference pattern appears.\n")

    print("Final Summary:")
    print("If D1 or D2 detect a photon, the result at D0 will not show an interference pattern.")
    print("If D3 or D4 detect a photon, the result at D0 will show an interference pattern.")
    print("\nThis corresponds to answer choice A.")


quantum_eraser_logic()