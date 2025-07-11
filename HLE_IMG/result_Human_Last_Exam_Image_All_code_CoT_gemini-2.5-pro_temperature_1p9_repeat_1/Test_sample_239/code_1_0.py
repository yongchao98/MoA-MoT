def analyze_quantum_eraser():
    """
    Analyzes the outcome of the delayed-choice quantum eraser experiment.

    This function determines whether an interference pattern will be observed at
    detector D0 based on which of the detectors D1, D2, D3, or D4 fires in
    coincidence.
    """

    # Rule: An interference pattern emerges if and only if which-path information is erased.
    rule = "Interference appears at D0 only if the 'which-path' information of its entangled partner is erased."

    # Detector Analysis
    # D1 and D2 provide which-path information because there is only one path to each.
    # D3 and D4 do not provide which-path information because two paths are recombined
    # at beam splitter BSc, making the original path indistinguishable.
    analysis = {
        'D1': {'path_info': 'Known', 'd0_outcome': 'No Interference Pattern'},
        'D2': {'path_info': 'Known', 'd0_outcome': 'No Interference Pattern'},
        'D3': {'path_info': 'Erased', 'd0_outcome': 'Interference Pattern'},
        'D4': {'path_info': 'Erased', 'd0_outcome': 'Interference Pattern'}
    }

    print("Analysis of the Delayed-Choice Quantum Eraser Experiment")
    print("-" * 60)
    print(f"Underlying Principle: {rule}")
    print("-" * 60)

    for detector, result in analysis.items():
        print(f"For a coincidence detection with detector {detector}:")
        print(f"  - Which-path information from the idler photon is '{result['path_info']}'.")
        # Here we 'output each number' by stating the result for each numbered detector.
        print(f"  - Result at D0: The subset of photons correlated with {detector} shows a '{result['d0_outcome']}'.")
        print("")

    print("Summary:")
    print("If D1 or D2 fires, the result at D0 will NOT show an interference pattern.")
    print("If D3 or D4 fires, the result at D0 WILL show an interference pattern.")
    print("\nThis logic corresponds to answer choice A.")


if __name__ == '__main__':
    analyze_quantum_eraser()