def solve_tnbc_question():
    """
    This function explains the rationale for the correct answer based on clinical trial data
    and then prints the final choice.
    """
    explanation = """
Based on major clinical trials, the addition of PD-1/PD-L1 inhibitors to chemotherapy for metastatic Triple Negative Breast Cancer (TNBC) has shown a significant benefit in a specific subgroup.

1.  **KEYNOTE-355 Trial:** This pivotal trial studied the PD-1 inhibitor pembrolizumab combined with chemotherapy. The results showed a statistically significant and clinically meaningful improvement in overall survival for patients whose tumors were PD-L1-positive (defined by a Combined Positive Score [CPS] ≥ 10).

2.  **IMpassion130 Trial:** Similarly, this trial studied the PD-L1 inhibitor atezolizumab with chemotherapy. It also demonstrated a significant overall survival benefit specifically in the PD-L1-positive patient population.

In both cases, the survival benefit was not statistically significant in the overall intention-to-treat (ITT) population, which includes all patients regardless of PD-L1 status. Furthermore, patients with PD-L1-negative tumors did not derive a significant survival benefit from this combination therapy.

Therefore, the population group in which PD-1 inhibitors, when added to chemotherapy, present a prolonged overall survival is the PD-L1-positive population.
"""
    print(explanation)
    final_answer = "B"
    print(f"The correct option is B.")
    print(f"<<<{final_answer}>>>")

solve_tnbc_question()