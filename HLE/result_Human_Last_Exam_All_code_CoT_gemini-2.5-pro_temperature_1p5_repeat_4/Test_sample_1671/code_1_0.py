import textwrap

def solve_medical_question():
    """
    This script analyzes the provided text to answer a specific question about NCCN guidelines.
    """
    # Step 1: Define the question and the provided answer choices.
    question = "Which of these has a category 1 recommendation from the NCCN for use in patients with HER2+ MBC whose disease progresses after treatment with trastuzumab and a taxane?"
    options = {
        "A": "Tucatinib trastuzumab and capecitabine",
        "B": "Trastuzumab emtansine",
        "C": "Fam-trastuzumab deruxtecan",
        "D": "Chemotherapy",
        "E": "None of the above"
    }

    # Step 2: Extract and analyze the relevant information from the text.
    # The text states the first-line treatment is THP (docetaxel, trastuzumab, pertuzumab).
    # Docetaxel is a taxane. So, the question is asking for the second-line treatment.
    relevant_text_nccn = "Current National Comprehensive Cancer Network (NCCN) guidelines with category 1 evidence recommend THP as the preferred first-line regimen for HER2+ MBC. Upon disease progression TDXd is recommended in the second-line setting."
    relevant_text_abbreviation = "Two HER2-directed ADCs ado-trastuzumab emtansine (T-DM1) and famtrastuzumab deruxtecan (TDXd) are approved for HER2+ MBC."

    print("Step-by-step analysis:")
    print("-----------------------")
    
    # Print the findings in a clear, wrapped format.
    print("1. Identifying the first-line treatment:")
    print(textwrap.fill("The question asks about treatment after progression on trastuzumab and a taxane. The text confirms the first-line NCCN recommended regimen is 'THP', which includes a taxane.", 80))
    
    print("\n2. Identifying the second-line treatment:")
    print(textwrap.fill(f"The text explicitly states the recommendation after first-line therapy: 'Upon disease progression TDXd is recommended in the second-line setting.'", 80))

    print("\n3. Identifying the full drug name for the abbreviation 'TDXd':")
    print(textwrap.fill(f"The text defines the abbreviation earlier: '...fam-trastuzumab deruxtecan (TDXd)...'", 80))

    # Step 3: Match the finding to the options.
    correct_treatment_name = "Fam-trastuzumab deruxtecan"
    final_answer_key = None
    for key, value in options.items():
        if value == correct_treatment_name:
            final_answer_key = key
            break
            
    print("\n4. Matching the drug to the provided options:")
    print(f"The drug 'Fam-trastuzumab deruxtecan' matches option {final_answer_key}.")
    print("-----------------------")

    # Step 4: Output the final answer in the required format.
    if final_answer_key:
        print("\nThe most likely diagnosis based on the NCCN guidelines cited in the text is:")
        print(f'<<<{final_answer_key}>>>')
    else:
        print("Could not determine the answer from the text.")

solve_medical_question()