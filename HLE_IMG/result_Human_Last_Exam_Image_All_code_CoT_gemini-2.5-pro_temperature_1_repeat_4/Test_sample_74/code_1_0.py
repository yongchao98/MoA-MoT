import numpy as np

def analyze_ngsc_data():
    """
    Analyzes NGSC data estimated from Figure 3b to evaluate answer choice L.
    """
    # Data estimated visually from Figure 3b, right panel.
    # Each participant has a list of NGSC values for 'psilocybin' (red dots)
    # and a single value for 'max_other' (the highest grey or blue dot).
    ngsc_data = {
        'P1': {'psilocybin': [0.75, 0.74, 0.74, 0.72], 'max_other': 0.72},
        'P3': {'psilocybin': [0.73, 0.73, 0.73, 0.71], 'max_other': 0.71},
        'P4': {'psilocybin': [0.82, 0.74, 0.74, 0.67], 'max_other': 0.70},
        'P5': {'psilocybin': [0.73, 0.72, 0.72, 0.68], 'max_other': 0.70},
        'P6': {'psilocybin': [0.74, 0.72], 'max_other': 0.68},
        'P7': {'psilocybin': [0.81, 0.80, 0.78, 0.75], 'max_other': 0.71}
    }

    print("Evaluating statement L: 'For all participants, mean NGSC in psilocybin condition was strictly higher than the maximum NGSC from any other condition.'")
    print("-" * 80)

    all_participants_hold = True
    for participant, data in ngsc_data.items():
        mean_psilocybin_ngsc = np.mean(data['psilocybin'])
        max_other_ngsc = data['max_other']
        
        # The comparison from statement L
        is_higher = mean_psilocybin_ngsc > max_other_ngsc
        
        print(f"Participant {participant}:")
        print(f"  Mean Psilocybin NGSC: {mean_psilocybin_ngsc:.4f}")
        print(f"  Maximum Other NGSC:   {max_other_ngsc:.4f}")
        # The statement is: mean(psilo) > max(other)
        print(f"  Is {mean_psilocybin_ngsc:.4f} > {max_other_ngsc:.4f}?  {is_higher}")
        print("-" * 20)
        
        if not is_higher:
            all_participants_hold = False

    print("\nFinal Conclusion:")
    if all_participants_hold:
        print("The statement holds true for all participants shown in the figure.")
        print("Therefore, answer choice L is directly supported by the data.")
    else:
        print("The statement does not hold true for all participants.")

analyze_ngsc_data()