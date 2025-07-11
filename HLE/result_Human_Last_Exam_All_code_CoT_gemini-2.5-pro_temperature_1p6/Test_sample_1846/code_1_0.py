import textwrap

def analyze_tnbc_treatment_outcomes():
    """
    Analyzes and explains the population group that benefits most from PD-1 inhibitors in TNBC.
    """
    
    explanation = """
    The question asks in which population group PD-1 inhibitors, when added to chemotherapy, show a prolonged overall survival for patients with Triple Negative Breast Cancer (TNBC).

    To answer this, we look at the results from major clinical trials that led to the approval and use of these drugs.

    1.  **PD-L1 as a Biomarker:** The protein PD-L1 (Programmed Death-Ligand 1) is a key predictive biomarker for the efficacy of PD-1/PD-L1 inhibitors. Its presence on tumor cells and immune cells suggests the cancer is using the PD-1/PD-L1 pathway to evade the immune system, making it a good target for these drugs.

    2.  **Clinical Trial Evidence (e.g., KEYNOTE-355):** The KEYNOTE-355 trial investigated the addition of pembrolizumab (a PD-1 inhibitor) to chemotherapy. The study found a statistically significant and clinically meaningful improvement in overall survival for patients whose tumors were PD-L1-positive (defined by a Combined Positive Score [CPS] ≥ 10).

    3.  **Outcome in Other Groups:**
        *   **Intention-to-treat (ITT) population:** This group includes all patients, regardless of PD-L1 status. The overall survival benefit in the full ITT population was not statistically significant.
        *   **PD-L1-negative population:** This group showed no significant survival benefit from the addition of the PD-1 inhibitor.

    **Conclusion:** The evidence strongly indicates that the prolonged overall survival benefit is primarily concentrated in the PD-L1-positive population.
    """

    answer_choice = "B"
    answer_text = "PD-L1-positive population"

    print("--- Medical Rationale ---")
    print(textwrap.dedent(explanation))
    print("--- Final Answer ---")
    print(f"The population group where PD-1 inhibitors show prolonged overall survival is the: ({answer_choice}) {answer_text}")

analyze_tnbc_treatment_outcomes()