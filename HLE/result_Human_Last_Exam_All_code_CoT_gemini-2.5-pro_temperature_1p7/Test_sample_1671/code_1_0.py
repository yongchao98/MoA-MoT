import textwrap

def find_treatment_recommendation():
    """
    Analyzes the provided medical text to find the NCCN recommended treatment
    for HER2+ MBC after progression on first-line therapy.
    """

    question = "Which of these has a category 1 recommendation from the NCCN for use in patients with HER2+ MBC whose disease progresses after treatment with trastuzumab and a taxane?"

    options = {
        'A': "Tucatinib trastuzumab and capecitabine",
        'B': "Trastuzumab emtansine",
        'C': "Fam-trastuzumab deruxtecan",
        'D': "Chemotherapy",
        'E': "None of the above"
    }

    # Information extracted from the provided text:
    first_line_therapy = "The text states that the NCCN recommended first-line therapy is THP (docetaxel [a taxane], trastuzumab, and pertuzumab)."
    second_line_therapy = "The text then states, 'Upon disease progression TDXd is recommended in the second-line setting.' and 'In 2022 TDXd replaced T-DM1 as second-line therapy...'"
    t_dxd_full_name = "TDXd is the abbreviation for fam-trastuzumab deruxtecan."
    conclusion = "Therefore, fam-trastuzumab deruxtecan is the recommended second-line therapy after progression on a regimen containing a taxane and trastuzumab."

    correct_option_key = 'C'
    
    # Print the analysis
    print("Analysis based on the provided text:")
    print("-" * 35)
    print(textwrap.fill(f"Question: {question}", width=80))
    print("\nRelevant findings from the text:")
    print(f"1. First-line treatment for HER2+ MBC: Trastuzumab and a taxane (as part of THP regimen).")
    print(f"2. After progression on first-line treatment, the recommended second-line treatment is TDXd (fam-trastuzumab deruxtecan).")
    
    print("\nConclusion:")
    print(f"The correct option is '{correct_option_key}', which is '{options[correct_option_key]}'.")


find_treatment_recommendation()