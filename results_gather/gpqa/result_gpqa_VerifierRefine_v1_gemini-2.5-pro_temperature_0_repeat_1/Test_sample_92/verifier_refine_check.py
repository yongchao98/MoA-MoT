import numpy as np
import math

def check_qpcr_correctness():
    """
    Checks the correctness of the LLM's answer by analyzing the qPCR data.
    """
    # Data from the problem description
    data = {
        100000: [33.9, 34.2, 34.5],
        10000: [30.6, 30.9, 31.2],
        1000: [27.3, 27.6, 27.9],
        100: [24, 24.3, 24.6],
        10: [20.7, 21, 21.3]
    }
    
    # The final answer provided by the LLM
    llm_answer = "A"

    # --- Step 1: Check Constraint A ---
    # "Ct values are not in agreement with the amount of target nucleic acid in samples"
    # Principle: Higher concentration should have a lower Ct value.
    concentrations = sorted(data.keys()) # Sort from low to high concentration
    avg_cts = [np.mean(data[c]) for c in concentrations]
    
    # For a correct experiment, as concentration increases, avg_ct should decrease.
    # The list `avg_cts` should be monotonically decreasing.
    is_relationship_inverted = False
    for i in range(len(avg_cts) - 1):
        # Check if Ct increases as concentration increases (the opposite of what's expected)
        if avg_cts[i] < avg_cts[i+1]:
            is_relationship_inverted = True
        else:
            # If it ever decreases, the trend is not consistently inverted.
            is_relationship_inverted = False
            break
            
    if not is_relationship_inverted:
        return "Incorrect. The analysis of constraint A is flawed. The provided data shows a direct relationship (higher concentration leads to higher Ct), which is the opposite of the expected inverse relationship in qPCR. The LLM correctly identified this as a major discrepancy, but the check failed to confirm it."

    # --- Step 2: Check Constraint C ---
    # "Ten-fold dilution is more than 3.3 cycles"
    # Let's check the difference between the average Ct of consecutive 10-fold dilutions.
    desc_concentrations = sorted(data.keys(), reverse=True)
    is_delta_ct_more_than_3_3 = False
    for i in range(len(desc_concentrations) - 1):
        ct1 = np.mean(data[desc_concentrations[i]])
        ct2 = np.mean(data[desc_concentrations[i+1]])
        delta_ct = ct1 - ct2
        if round(delta_ct, 2) > 3.3:
            is_delta_ct_more_than_3_3 = True
            break
            
    if is_delta_ct_more_than_3_3:
        return "Incorrect. The analysis of constraint C is flawed. The statement is 'Ten-fold dilution is more than 3.3 cycles'. The data shows the difference is exactly 3.3 cycles for each dilution step. Therefore, the statement is false, which the LLM correctly identified."

    # --- Step 3: Check Constraint D ---
    # "The deviation is more than 0.3 between technical replicates"
    # Let's check the range (max - min) for each triplicate.
    is_deviation_high = True
    for conc, cts in data.items():
        deviation = max(cts) - min(cts)
        # The statement is true if the deviation is > 0.3. We check if this holds for all.
        if not (deviation > 0.3):
            is_deviation_high = False
            break
        # Let's also check if the LLM's calculation of 0.6 is correct
        if round(deviation, 2) != 0.6:
             return f"Incorrect. The LLM's calculation for deviation at {conc} is wrong. The code calculated {round(deviation,2)} but the LLM stated 0.6."

    if not is_deviation_high:
        return "Incorrect. The analysis of constraint D is flawed. The statement is 'The deviation is more than 0.3'. The data shows the deviation for all replicates is 0.6, which is indeed > 0.3. Therefore, the statement is true, which the LLM correctly identified."

    # --- Final Conclusion ---
    # The code confirms the following:
    # 1. Constraint A is TRUE: The relationship between concentration and Ct is inverted.
    # 2. Constraint C is FALSE: The delta Ct is exactly 3.3, not more.
    # 3. Constraint D is TRUE: The deviation within replicates is 0.6, which is > 0.3.
    #
    # The LLM correctly identifies that both A and D are true statements about the data.
    # It then correctly reasons that A represents a fundamental failure of the experiment's validity (likely due to mislabeling),
    # making it the most significant discrepancy that invalidates the entire standard curve.
    # The choice of A as the best explanation is therefore correct.
    
    if llm_answer == "A":
        return "Correct"
    else:
        return f"Incorrect. The correct answer is A. The most fundamental error is that the Ct values are inverted relative to the concentration, making the entire standard curve invalid. The LLM chose {llm_answer}."

# Run the check
result = check_qpcr_correctness()
print(result)