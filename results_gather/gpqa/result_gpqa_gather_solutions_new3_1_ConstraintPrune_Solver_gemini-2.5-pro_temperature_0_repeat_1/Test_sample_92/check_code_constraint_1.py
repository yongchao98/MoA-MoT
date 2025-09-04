import numpy as np

def check_qpcr_answer():
    """
    Checks the correctness of the provided answer for the qPCR analysis question.
    """
    # --- 1. Data and Question Setup ---
    
    # Data from the problem description
    data = {
        100000: [33.9, 34.2, 34.5],
        10000: [30.6, 30.9, 31.2],
        1000: [27.3, 27.6, 27.9],
        100: [24, 24.3, 24.6],
        10: [20.7, 21, 21.3]
    }

    # The final answer provided by the LLM to be checked.
    # The prompt's final analysis points to 'B'.
    llm_answer = "B"

    # Options from the question
    options = {
        "A": "Ten-fold dilution is more than 3.3 cycles",
        "B": "Ct values are not in agreement with the amount of target nucleic acid in samples",
        "C": "The deviation is more than 0.3 between technical replicates",
        "D": "qPCR cannot be used for the quantification of nucleic acid in samples"
    }

    # --- 2. Constraint and Option Verification ---
    
    analysis_results = {}

    # Check Option A: "Ten-fold dilution is more than 3.3 cycles"
    desc_concentrations = sorted(data.keys(), reverse=True)
    avg_cts = [np.mean(data[c]) for c in desc_concentrations]
    diffs = [avg_cts[i] - avg_cts[i+1] for i in range(len(avg_cts)-1)]
    # The statement is that the difference is MORE than 3.3. Let's check if any are.
    is_A_true = any(d > 3.3 for d in diffs)
    analysis_results['A'] = {
        "is_true": is_A_true,
        "reason": f"The differences between 10-fold dilutions are {[round(d, 2) for d in diffs]}. All are exactly 3.3, not more. So, statement A is false."
    }

    # Check Option B: "Ct values are not in agreement with the amount of target nucleic acid in samples"
    # This implies the fundamental inverse relationship is broken.
    # A correct relationship means as concentration increases, Ct decreases.
    asc_concentrations = sorted(data.keys())
    asc_avg_cts = [np.mean(data[c]) for c in asc_concentrations]
    # Check if Ct values are INCREASING with concentration (a direct relationship)
    is_direct_relationship = all(asc_avg_cts[i] < asc_avg_cts[i+1] for i in range(len(asc_avg_cts)-1))
    is_B_true = is_direct_relationship
    analysis_results['B'] = {
        "is_true": is_B_true,
        "reason": f"The relationship between concentration and Ct value is direct (higher concentration gives higher Ct), which contradicts the fundamental inverse principle of qPCR. So, statement B is true."
    }

    # Check Option C: "The deviation is more than 0.3 between technical replicates"
    all_deviations_over_0_3 = True
    deviation_reasons = []
    for conc, cts in data.items():
        deviation = max(cts) - min(cts)
        deviation_reasons.append(f"For {conc} copies, deviation is {deviation:.1f}.")
        if deviation <= 0.3:
            all_deviations_over_0_3 = False
    is_C_true = all_deviations_over_0_3
    analysis_results['C'] = {
        "is_true": is_C_true,
        "reason": f"The deviation (max-min) for all technical replicates is > 0.3. {' '.join(deviation_reasons)} So, statement C is true."
    }

    # Check Option D: "qPCR cannot be used for the quantification of nucleic acid in samples"
    # This is a general knowledge check.
    is_D_true = False
    analysis_results['D'] = {
        "is_true": is_D_true,
        "reason": "qPCR is a gold-standard, scientifically validated method for nucleic acid quantification. So, statement D is false."
    }

    # --- 3. Final Decision Logic ---
    
    # Identify all true statements about the data
    true_options = [opt for opt, result in analysis_results.items() if result["is_true"]]

    if not true_options:
        return "Error: No option correctly describes a discrepancy in the data."

    # Determine the BEST explanation. A fundamental error (B) is more significant than a precision error (C).
    best_explanation = None
    if 'B' in true_options:
        best_explanation = 'B'
    elif 'C' in true_options:
        best_explanation = 'C'
    elif 'A' in true_options:
        best_explanation = 'A'
    else:
        # This case should not be reached if there are true options.
        best_explanation = true_options[0]

    # --- 4. Compare with LLM's Answer and Return Result ---
    
    if llm_answer == best_explanation:
        return "Correct"
    else:
        reasoning = f"The provided answer '{llm_answer}' is incorrect.\n"
        reasoning += f"My analysis determined the best explanation is '{best_explanation}'.\n"
        reasoning += "Here is the breakdown:\n"
        reasoning += f" - Option A is false. {analysis_results['A']['reason']}\n"
        reasoning += f" - Option D is false. {analysis_results['D']['reason']}\n"
        reasoning += f" - Both B and C are factually true statements about the data:\n"
        reasoning += f"   - B: {analysis_results['B']['reason']}\n"
        reasoning += f"   - C: {analysis_results['C']['reason']}\n"
        reasoning += "However, the question asks for the best explanation. The inverted relationship described in B is a fundamental, systematic error that invalidates the entire experiment's conclusion. The poor precision described in C is a secondary issue. Therefore, B is the most significant discrepancy."
        return reasoning

# Execute the check and print the result
print(check_qpcr_answer())