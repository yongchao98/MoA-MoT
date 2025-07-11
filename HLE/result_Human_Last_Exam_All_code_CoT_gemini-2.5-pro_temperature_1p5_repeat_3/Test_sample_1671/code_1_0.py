import textwrap

def find_treatment_recommendation():
    """
    Analyzes the provided text to determine the NCCN-recommended second-line treatment
    for HER2+ MBC after progression on trastuzumab and a taxane.
    """
    
    # The question asks for the recommended second-line therapy for HER2+ MBC.
    question_context = "NCCN Category 1 recommendation for use in patients with HER2+ MBC whose disease progresses after treatment with trastuzumab and a taxane."

    # Analysis based on the provided text:
    # Step 1: Identify the first-line therapy to establish the starting point.
    first_line_therapy_quote = "Current National Comprehensive Cancer Network (NCCN) guidelines with category 1 evidence recommend THP as the preferred first-line regimen for HER2+ MBC."
    # Note: THP is docetaxel (a taxane), trastuzumab, and pertuzumab. So the question is asking for the next step.

    # Step 2: Identify the recommended therapy upon progression (i.e., second-line therapy).
    second_line_therapy_quote = "Upon disease progression, TDXd is recommended in the second-line setting."
    # Note: TDXd stands for fam-trastuzumab deruxtecan. This corresponds to option C.

    # Step 3: Check the status of other options mentioned in the text.
    t_dm1_status_quote = "In 2022, TDXd replaced T-DM1 as second-line therapy..."
    # Note: This explicitly disqualifies T-DM1 (Trastuzumab emtansine, option B) as the current standard.
    
    tucatinib_regimen_status_quote = "...NCCN guidelines recommend tucatinib, trastuzumab, and capecitabine as the preferred regimen [for] third-line therapy."
    # Note: This places the tucatinib regimen (option A) as a third-line option, not second-line.
    
    # Conclusion: The text clearly indicates fam-trastuzumab deruxtecan is the recommended second-line therapy.
    
    print("Based on the analysis of the provided text:")
    print("-" * 30)
    print(f"The question asks for the standard second-line therapy for HER2+ MBC after progression on a regimen like trastuzumab and a taxane.")
    print("\n1. The text establishes the standard first-line therapy: '...NCCN guidelines with category 1 evidence recommend THP [trastuzumab, pertuzumab, and a taxane] as the preferred first-line regimen...'")
    print("\n2. The text then specifies the next step: 'Upon disease progression, TDXd [fam-trastuzumab deruxtecan] is recommended in the second-line setting.'")
    print("\n3. This makes Fam-trastuzumab deruxtecan the correct answer for the second-line setting described in the question.")
    print("\n4. Other options are incorrect as per the text: Trastuzumab emtansine was replaced by TDXd, and the Tucatinib regimen is recommended for third-line therapy.")
    print("-" * 30)
    print("Therefore, the correct answer is C.")
    
    final_answer = "<<<C>>>"
    print(f"\nFinal Answer: {final_answer}")

find_treatment_recommendation()