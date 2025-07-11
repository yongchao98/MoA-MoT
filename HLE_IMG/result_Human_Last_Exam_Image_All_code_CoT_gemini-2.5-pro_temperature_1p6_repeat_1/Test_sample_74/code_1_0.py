import numpy as np
from scipy import stats

def analyze_ngsc_data():
    """
    Analyzes estimated NGSC data from Figure 3b to verify answer choice G.
    """
    # Visually estimated Whole-brain NGSC values from Figure 3b (right panel).
    # 'no_drug_ngsc' corresponds to the grey dots (baseline/"pre-psilocybin").
    no_drug_ngsc = [
        # P1
        0.71, 0.69, 0.68, 0.70, 0.67, 0.66,
        # P3
        0.69, 0.67, 0.71, 0.66, 0.68, 0.70,
        # P4
        0.66, 0.65, 0.68, 0.67, 0.67,
        # P5
        0.66, 0.63, 0.62, 0.65, 0.67, 0.64,
        # P6 (data from both blocks)
        0.70, 0.69, 0.67, 0.63, 0.67, 0.66, 0.68, 0.67, 0.66, 0.66,
        # P7
        0.69, 0.68, 0.65, 0.67, 0.70, 0.66
    ]

    # 'psilocybin_ngsc' corresponds to the red dots ("post-psilocybin").
    psilocybin_ngsc = [
        # P1
        0.73, 0.76, 0.70, 0.68,
        # P3
        0.73, 0.72, 0.75, 0.71,
        # P4
        0.71, 0.68, 0.72, 0.70,
        # P5
        0.73, 0.70, 0.74, 0.69,
        # P6
        0.72, 0.74, 0.71,
        # P7
        0.79, 0.81, 0.77, 0.76
    ]

    # --- Components for the t-test equation ---
    # t = (mean1 - mean2) / sqrt(var1/n1 + var2/n2)
    
    # Group 1: Psilocybin
    mean1 = np.mean(psilocybin_ngsc)
    var1 = np.var(psilocybin_ngsc, ddof=1) # ddof=1 for sample variance
    n1 = len(psilocybin_ngsc)

    # Group 2: No Drug
    mean2 = np.mean(no_drug_ngsc)
    var2 = np.var(no_drug_ngsc, ddof=1)
    n2 = len(no_drug_ngsc)

    print("--- T-Test Equation Components ---")
    print(f"Psilocybin Group Mean (mean1): {mean1:.4f}")
    print(f"Psilocybin Group Sample Size (n1): {n1}")
    print(f"Psilocybin Group Variance (var1): {var1:.6f}\n")

    print(f"No Drug Group Mean (mean2): {mean2:.4f}")
    print(f"No Drug Group Sample Size (n2): {n2}")
    print(f"No Drug Group Variance (var2): {var2:.6f}\n")
    
    # --- Perform Independent Samples T-test ---
    t_statistic, p_value = stats.ttest_ind(psilocybin_ngsc, no_drug_ngsc, equal_var=False)

    print("--- T-Test Results ---")
    print(f"T-statistic: {t_statistic:.4f}")
    print(f"P-value: {p_value}")

    print("\n--- Conclusion ---")
    if p_value < 0.05:
        print("The p-value is much less than 0.05, indicating a statistically significant difference between the two groups.")
        print("The positive t-statistic shows the mean of the psilocybin group is significantly higher than the no-drug group.")
        print("This computational result confirms the visual evidence from the graph and directly supports answer choice G.")
    else:
        print("The result is not statistically significant.")

analyze_ngsc_data()
<<<G>>>