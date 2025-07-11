def explain_quantum_eraser():
    """
    Explains the expected results for the delayed-choice quantum eraser experiment.
    """
    
    # Principle of Quantum Interference
    interference_principle = "Interference occurs if and only if 'which-path' information is unknowable."
    no_interference_principle = "If 'which-path' information is known, interference is destroyed."

    print("Analyzing the Delayed-Choice Quantum Eraser Experiment:\n")
    print(f"The core principle is: {interference_principle}\n")

    # Analysis for Detectors D1 and D2
    print("--- Case 1: Detection at D1 or D2 ---")
    print("A photon detected at D1 unambiguously came from one slit.")
    print("A photon detected at D2 unambiguously came from the other slit.")
    print("Therefore, a detection at D1 or D2 provides 'which-path' information.")
    print(f"Result based on principle: {no_interference_principle}")
    print("Conclusion: For D0 data correlated with D1 or D2, NO interference pattern will be shown.\n")

    # Analysis for Detectors D3 and D4
    print("--- Case 2: Detection at D3 or D4 ---")
    print("Paths from both slits are recombined at beam-splitter BSc before reaching D3 and D4.")
    print("A photon at D3 or D4 could have originated from either slit.")
    print("Therefore, 'which-path' information is erased.")
    print(f"Result based on principle: {interference_principle}")
    print("Conclusion: For D0 data correlated with D3 or D4, an interference pattern WILL be shown.\n")

    # Final Summary
    print("--- Final Summary ---")
    final_answer_text = "If D3 or D4 detect a photon, the result at D0 will show an interference pattern. If D1 or D2 detect a photon, the result at D0 will not show an interference pattern."
    print(final_answer_text)

explain_quantum_eraser()