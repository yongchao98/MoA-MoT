def find_tnbc_treatment_population():
    """
    This script explains the population group in which PD-1 inhibitors
    showed prolonged overall survival for Triple Negative Breast Cancer (TNBC)
    based on major clinical trial results.
    """

    explanation = """
Clinical trials, such as KEYNOTE-355, have been pivotal in determining the efficacy of PD-1 inhibitors (like pembrolizumab) in combination with chemotherapy for treating metastatic Triple Negative Breast Cancer (TNBC).

The findings of this study can be summarized as follows:
- Intention-to-treat (ITT) population: This group includes all randomized patients, regardless of their biomarker status. In the KEYNOTE-355 trial, the improvement in overall survival for the ITT population was not statistically significant.
- PD-L1-negative population: Patients in this group did not show a significant survival benefit from the addition of a PD-1 inhibitor.
- PD-L1-positive population: A statistically significant and clinically meaningful improvement in overall survival was observed specifically in patients with PD-L1-positive tumors (especially those with a high expression, CPS ≥ 10).

Therefore, the benefit of prolonged overall survival is primarily seen in the PD-L1-positive population, making PD-L1 expression a critical predictive biomarker for this therapy.
"""

    answer_choice = "B"
    answer_text = "PD-L1-positive population"

    print(explanation)
    print(f"Conclusion: The correct choice is '{answer_choice}', which corresponds to the '{answer_text}'.")

find_tnbc_treatment_population()