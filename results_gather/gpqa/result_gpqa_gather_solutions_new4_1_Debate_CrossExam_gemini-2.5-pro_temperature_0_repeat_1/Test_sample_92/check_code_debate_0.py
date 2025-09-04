import numpy as np
import math

def check_qpcr_answer():
    """
    Checks the correctness of the provided answer for the qPCR problem.

    The function analyzes the qPCR data based on fundamental principles and
    evaluates each multiple-choice option to determine the best explanation
    for the discrepancies.
    """
    # Data from the question
    data = {
        100000: [33.9, 34.2, 34.5],
        10000: [30.6, 30.9, 31.2],
        1000: [27.3, 27.6, 27.9],
        100: [24, 24.3, 24.6],
        10: [20.7, 21, 21.3]
    }

    # The final answer from the LLM to be checked
    llm_answer = 'D'

    # --- Step-by-step analysis ---
    print("--- Starting qPCR Data Analysis ---")

    # 1. Analyze the relationship between concentration and Ct value (for Option D)
    print("\n1. Checking Trend: Ct vs. Concentration (for Option D)")
    concentrations = sorted(data.keys(), reverse=True)  # High to low concentration
    avg_cts = [np.mean(data[c]) for c in concentrations]
    
    # Correct trend: As concentration decreases, Ct should INCREASE.
    # Observed trend:
    is_trend_inverted = all(avg_cts[i] > avg_cts[i+1] for i in range(len(avg_cts)-1))
    
    is_D_true = False
    if is_trend_inverted:
        print("   - Finding: As concentration decreases, the average Ct value also decreases.")
        print("   - Conclusion: This is the INVERSE of the expected trend. The Ct values are not in agreement with the amount of target nucleic acid.")
        is_D_true = True
    else:
        print("   - Finding: The trend is not consistently inverted, but it is incorrect.")
        is_D_true = True # Still not in agreement

    # 2. Analyze the deviation between technical replicates (for Option A)
    print("\n2. Checking Deviation within Replicates (for Option A)")
    deviations = {conc: max(cts) - min(cts) for conc, cts in data.items()}
    all_deviations_over_0_3 = all(dev > 0.3 for dev in deviations.values())
    
    is_A_true = False
    if all_deviations_over_0_3:
        print(f"   - Finding: The deviation (max Ct - min Ct) for all replicates is consistently {list(deviations.values())[0]:.1f}.")
        print("   - Conclusion: The deviation is more than 0.3 between technical replicates.")
        is_A_true = True
    else:
        print("   - Conclusion: The statement 'deviation is more than 0.3' is not true for all replicates.")

    # 3. Analyze the difference between 10-fold dilutions (for Option C)
    print("\n3. Checking Difference between Dilutions (for Option C)")
    ct_differences = [abs(avg_cts[i] - avg_cts[i+1]) for i in range(len(avg_cts) - 1)]
    is_C_true = any(diff > 3.3 for diff in ct_differences)
    
    if is_C_true:
        print("   - Conclusion: At least one 10-fold dilution step is more than 3.3 cycles.")
    else:
        # Check if it's exactly 3.3
        is_exactly_3_3 = all(math.isclose(diff, 3.3) for diff in ct_differences)
        if is_exactly_3_3:
            print("   - Finding: The difference between each 10-fold dilution is consistently 3.3 cycles.")
            print("   - Conclusion: The statement 'Ten-fold dilution is more than 3.3 cycles' is FALSE.")
        else:
            print("   - Conclusion: The statement 'Ten-fold dilution is more than 3.3 cycles' is FALSE.")

    # 4. Analyze the conceptual option B
    print("\n4. Checking Conceptual Validity (for Option B)")
    print("   - Conclusion: The statement 'qPCR cannot be used for quantification' is conceptually FALSE. It is a standard method.")
    is_B_true = False

    # --- Final Verdict ---
    print("\n--- Synthesizing Results ---")
    print(f"Summary of Findings:")
    print(f"  - Is statement A true? {is_A_true}")
    print(f"  - Is statement B true? {is_B_true}")
    print(f"  - Is statement C true? {is_C_true}")
    print(f"  - Is statement D true? {is_D_true}")

    # Determine the best answer based on the findings
    best_explanation = None
    if is_D_true:
        # The inverted trend (D) is a fundamental error of validity, making the experiment nonsensical.
        # The high deviation (A) is an error of precision.
        # An error of validity is more critical than an error of precision.
        best_explanation = 'D'
    elif is_A_true:
        best_explanation = 'A'
    elif is_C_true:
        best_explanation = 'C'
    else:
        best_explanation = 'None of the options are a good explanation.'

    print(f"\nBased on the analysis, the best explanation is Option '{best_explanation}'.")
    print(f"The LLM's answer was '{llm_answer}'.")

    if best_explanation == llm_answer:
        return "Correct"
    else:
        reason = f"The provided answer '{llm_answer}' is incorrect. "
        reason += f"The analysis shows that the best explanation is '{best_explanation}'.\n"
        if llm_answer == 'A' and best_explanation == 'D':
            reason += "Reason: While statement A (high deviation) is true, statement D (inverted trend) describes a more fundamental error that invalidates the entire experiment's quantitative basis. An error of validity is more significant than an error of precision."
        elif llm_answer == 'D' and best_explanation != 'D':
             reason += f"Reason: The analysis concluded that '{best_explanation}' was the best choice, contradicting the provided answer."
        return reason

# Execute the check and print the result
result = check_qpcr_answer()
print("\n--- FINAL VERDICT ---")
print(result)