import textwrap

def explain_tnbc_immunotherapy():
    """
    Explains the basis for choosing the correct population group for PD-1 inhibitor efficacy in TNBC.
    """
    explanation = """
    The question asks in which population group PD-1 inhibitors, when added to chemotherapy, lead to prolonged overall survival in patients with Triple Negative Breast Cancer (TNBC). The answer is based on data from pivotal clinical trials.

    1.  **Biomarker Importance**: The expression of PD-L1 (Programmed death-ligand 1) on tumor cells and immune cells is a key biomarker used to predict which patients are more likely to respond to immunotherapy drugs like PD-1 inhibitors.

    2.  **Key Clinical Trial (KEYNOTE-355)**: This major trial studied the effect of adding pembrolizumab (a PD-1 inhibitor) to chemotherapy for metastatic TNBC.

    3.  **Trial Results**: The results showed that the greatest benefit was concentrated in patients whose tumors expressed PD-L1. A statistically significant and clinically meaningful improvement in overall survival (OS) was observed in the PD-L1-positive population (specifically those with a Combined Positive Score [CPS] ≥ 10).

    4.  **Other Populations**: While some benefit for other metrics (like progression-free survival) was seen in the broader intention-to-treat population, the significant advantage in overall survival—the gold standard for cancer treatment efficacy—was clearly established in the PD-L1-positive subgroup. The PD-L1-negative population did not derive a significant survival benefit from the addition of the PD-1 inhibitor.

    5.  **Conclusion**: Therefore, the PD-L1-positive population is the group where PD-1 inhibitors have been shown to prolong overall survival compared to chemotherapy alone.
    """
    print(textwrap.dedent(explanation).strip())

explain_tnbc_immunotherapy()
print("\n<<<B>>>")