import numpy as np

def analyze_ngsc_data():
    """
    Analyzes estimated NGSC data from Siegel et al. 2024, Figure 3b,
    to verify the effect of psilocybin.
    """
    # Data is visually estimated from Figure 3b (right panel).
    # We will compare the mean of psilocybin scans vs. the combined mean of
    # control scans (MTP and No Drug).
    ngsc_data = {
        'P1': {
            'control': [0.71, 0.69, 0.68, 0.65, 0.72],
            'psilocybin': [0.73, 0.74, 0.75, 0.76]
        },
        'P3': {
            'control': [0.63, 0.66, 0.67, 0.70, 0.71],
            'psilocybin': [0.71, 0.73, 0.75]
        },
        'P4': {
            'control': [0.63, 0.65, 0.67, 0.70, 0.71],
            'psilocybin': [0.67, 0.74, 0.82, 0.75]
        },
        'P5': {
            'control': [0.63, 0.66, 0.67, 0.70],
            'psilocybin': [0.69, 0.71, 0.72, 0.75]
        },
        'P6': {
            'control': [0.62, 0.65, 0.64, 0.67, 0.68, 0.69],
            'psilocybin': [0.715, 0.72, 0.74]
        },
        'P7': {
            'control': [0.63, 0.66, 0.68, 0.71],
            'psilocybin': [0.74, 0.78, 0.79, 0.81]
        }
    }

    print("Analyzing NGSC data based on Figure 3b from Siegel et al. 2024.")
    print("Hypothesis (Choice G): Post-psilocybin whole-brain NGSC shows a significant increase compared to pre-psilocybin (control) NGSC.\n")

    all_participants_show_increase = True
    for participant, conditions in ngsc_data.items():
        mean_control = np.mean(conditions['control'])
        mean_psilocybin = np.mean(conditions['psilocybin'])

        # The final equation here is a simple comparison of the means.
        # We output each number in this comparison.
        print(f"Participant {participant}:")
        print(f"  - Mean Control NGSC: {mean_control:.4f}")
        print(f"  - Mean Psilocybin NGSC: {mean_psilocybin:.4f}")

        if mean_psilocybin > mean_control:
            print(f"  - Result: The mean Psilocybin NGSC ({mean_psilocybin:.4f}) is greater than the mean Control NGSC ({mean_control:.4f}).")
        else:
            print(f"  - Result: The mean Psilocybin NGSC ({mean_psilocybin:.4f}) is NOT greater than the mean Control NGSC ({mean_control:.4f}).")
            all_participants_show_increase = False
        print("-" * 20)

    print("\nConclusion:")
    if all_participants_show_increase:
        print("For all participants analyzed, the mean NGSC under the psilocybin condition is higher than under control conditions.")
        print("This computational analysis of estimated data supports answer choice G.")
    else:
        print("The analysis does not consistently support answer choice G for all participants.")

if __name__ == '__main__':
    analyze_ngsc_data()