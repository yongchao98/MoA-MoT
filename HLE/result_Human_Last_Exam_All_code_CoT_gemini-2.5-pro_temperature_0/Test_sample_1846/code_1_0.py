def solve_tnbc_question():
    """
    This function explains the clinical trial evidence regarding PD-1 inhibitors in
    Triple Negative Breast Cancer (TNBC) and provides the correct answer.
    """
    explanation = """
Clinical trials have shown that the benefit of adding PD-1/PD-L1 inhibitors to chemotherapy for metastatic Triple Negative Breast Cancer (TNBC) is most significant in a specific subgroup of patients.

1.  **KEYNOTE-355 Trial:** This trial studied pembrolizumab plus chemotherapy. It demonstrated a statistically significant and clinically meaningful improvement in both progression-free survival and overall survival for patients with PD-L1-positive tumors (specifically, those with a Combined Positive Score [CPS] ≥ 10). The benefit for overall survival was not statistically significant in the broader intention-to-treat (ITT) population.

2.  **IMpassion130 Trial:** This trial studied atezolizumab plus chemotherapy. It showed a significant overall survival benefit in the pre-specified analysis for the PD-L1-positive population. However, this benefit did not reach statistical significance in the final analysis for the overall intention-to-treat (ITT) population.

**Conclusion:**
Based on the results of these major trials, the addition of a PD-1/PD-L1 inhibitor to chemotherapy provides a prolonged overall survival benefit specifically for the **PD-L1-positive population** when compared to chemotherapy alone.
"""
    print(explanation)
    final_answer = "B"
    print(f"<<<{final_answer}>>>")

solve_tnbc_question()