import numpy as np
from scipy import stats

def analyze_ngsc_data():
    """
    Analyzes estimated NGSC data from Figure 3b to verify answer choice G.
    """
    # Data estimated visually from Figure 3b, right panel
    # (P2 is missing from the original figure's x-axis)
    ngsc_data = {
        'P1': {
            'Psilocybin': [0.73, 0.74, 0.67],
            'MTP': [0.67, 0.69, 0.70],
            'No drug': [0.65, 0.69, 0.71, 0.71, 0.72]
        },
        'P3': {
            'Psilocybin': [0.69, 0.72, 0.74],
            'MTP': [0.67, 0.68, 0.69],
            'No drug': [0.62, 0.64, 0.67, 0.68, 0.70, 0.71]
        },
        'P4': {
            'Psilocybin': [0.70, 0.74, 0.82],
            'MTP': [0.66, 0.67, 0.69, 0.71],
            'No drug': [0.65, 0.67, 0.68, 0.68, 0.69, 0.70, 0.72]
        },
        'P5': {
            'Psilocybin': [0.65, 0.70, 0.73],
            'MTP': [0.67, 0.68, 0.68, 0.69],
            'No drug': [0.64, 0.65, 0.66, 0.68, 0.69, 0.70, 0.71]
        },
        'P6': {
            'Psilocybin': [0.71, 0.73],
            'MTP': [0.64, 0.65, 0.66, 0.67],
            'No drug': [0.63, 0.66, 0.67, 0.67, 0.68, 0.69, 0.70]
        },
        'P7': {
            'Psilocybin': [0.76, 0.78, 0.79, 0.80, 0.81],
            'MTP': [0.68, 0.69, 0.70, 0.72],
            'No drug': [0.66, 0.68, 0.69, 0.69, 0.70, 0.71]
        }
    }
    
    print("Analysis of NGSC Increase for each participant:\n")
    all_participants_pre = []
    all_participants_psilo = []

    for participant, data in ngsc_data.items():
        pre_psilocybin_scans = data['MTP'] + data['No drug']
        psilocybin_scans = data['Psilocybin']
        
        all_participants_pre.extend(pre_psilocybin_scans)
        all_participants_psilo.extend(psilocybin_scans)
        
        # Calculate means
        mean_pre = np.mean(pre_psilocybin_scans)
        mean_psilo = np.mean(psilocybin_scans)
        
        # Perform Welch's t-test (assumes unequal variances)
        t_stat, p_value = stats.ttest_ind(psilocybin_scans, pre_psilocybin_scans, equal_var=False)
        
        print(f"--- Participant {participant} ---")
        print(f"Mean Pre-Psilocybin NGSC: {mean_pre:.4f}")
        print(f"Mean Psilocybin NGSC:     {mean_psilo:.4f}")
        print(f"T-test result: p-value = {p_value:.4f}")
        if p_value < 0.05 and mean_psilo > mean_pre:
            print("Conclusion: Significant increase in NGSC observed.")
        else:
            print("Conclusion: Increase not statistically significant at p < 0.05 (or is a decrease).")
        print("-" * 25 + "\n")

    # Overall analysis
    mean_all_pre = np.mean(all_participants_pre)
    mean_all_psilo = np.mean(all_participants_psilo)
    t_stat_all, p_value_all = stats.ttest_ind(all_participants_psilo, all_participants_pre, equal_var=False)

    print("--- Overall Analysis Across All Participants ---")
    print(f"Mean Pre-Psilocybin NGSC: {mean_all_pre:.4f}")
    print(f"Mean Psilocybin NGSC:     {mean_all_psilo:.4f}")
    print(f"T-test result: p-value = {p_value_all:.4e}") # Use scientific notation for very small p-value
    if p_value_all < 0.05 and mean_all_psilo > mean_all_pre:
        print("Conclusion: A highly significant increase in NGSC is observed overall.")
    else:
        print("Conclusion: The overall increase is not statistically significant.")

if __name__ == '__main__':
    analyze_ngsc_data()