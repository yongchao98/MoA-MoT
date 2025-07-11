def explain_tnbc_treatment_population():
    """
    This function explains which patient population with Triple Negative Breast Cancer (TNBC)
    benefits from PD-1 inhibitors based on clinical trial data.
    """
    explanation = """
Based on data from major clinical trials, the efficacy of PD-1 inhibitors in combination with chemotherapy for Triple Negative Breast Cancer (TNBC) is closely linked to the expression of a specific biomarker.

1.  **Key Biomarker:** The crucial biomarker is Programmed Death-Ligand 1 (PD-L1). Patients are often categorized as 'PD-L1-positive' or 'PD-L1-negative' based on the level of this protein in and around the tumor.

2.  **Clinical Trial Evidence:** The KEYNOTE-355 trial, a landmark study, evaluated the PD-1 inhibitor pembrolizumab plus chemotherapy. It found a statistically significant and clinically meaningful improvement in overall survival.

3.  **Benefiting Population:** This survival benefit was specifically observed in the **PD-L1-positive population**. In contrast, the overall survival benefit was not statistically significant in the intention-to-treat (ITT) population, which includes both PD-L1-positive and PD-L1-negative patients. The lack of benefit in the PD-L1-negative subgroup diluted the positive effect seen in the overall group.

Therefore, the PD-L1-positive population is the group where PD-1 inhibitor treatment presents a prolonged overall survival compared to chemotherapy alone.
"""
    print(explanation)

explain_tnbc_treatment_population()