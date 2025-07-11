def find_tnbc_treatment_answer():
    """
    This function provides the answer and explanation regarding the efficacy
    of PD-1 inhibitors in Triple Negative Breast Cancer (TNBC) treatment.
    """
    
    question = "In comparison to chemotherapy alone, in which population group the PD-1 inhibitors treatment presents a prolonged overall survival for Triple Negative Breast Cancer (TNBC)?"
    
    options = {
        'A': 'Intention-to-treat population',
        'B': 'PD-L1-positive population',
        'C': 'Both intention-to-treat population and PD-L1-positive population',
        'D': 'PD-L1-negative population',
        'E': 'None of above'
    }

    correct_answer = 'B'
    
    explanation = (
        "Based on landmark clinical trials such as KEYNOTE-355, the addition of a PD-1 inhibitor (pembrolizumab) "
        "to chemotherapy demonstrated a statistically significant and clinically meaningful improvement "
        "in overall survival (OS) compared to chemotherapy alone in patients with previously untreated "
        "metastatic Triple Negative Breast Cancer whose tumors expressed PD-L1 (with a Combined Positive Score [CPS] of 10 or more).\n\n"
        "This survival benefit was not statistically significant in the intention-to-treat (ITT) population as a whole, "
        "which includes patients regardless of their PD-L1 status. Therefore, the benefit is primarily confined to the "
        "PD-L1-positive population."
    )

    print(f"Question: {question}\n")
    print("Answer Choices:")
    for key, value in options.items():
        print(f"  {key}. {value}")
    
    print("\n" + "="*30)
    print("   Final Answer and Justification   ")
    print("="*30)
    print(f"\nThe correct answer is: {correct_answer} ({options[correct_answer]})")
    print("\nExplanation:")
    print(explanation)

# Execute the function to display the answer
find_tnbc_treatment_answer()