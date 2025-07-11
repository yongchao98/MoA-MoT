import textwrap

def solve_medical_question():
    """
    Analyzes the provided text to answer a multiple-choice question about cancer treatment.
    """
    question = "Which of these has a category 1 recommendation from the NCCN for use in patients with HER2+ MBC whose disease progresses after treatment with trastuzumab and a taxane?"
    
    options = {
        'A': 'Tucatinib trastuzumab and capecitabine',
        'B': 'Trastuzumab emtansine',
        'C': 'Fam-trastuzumab deruxtecan',
        'D': 'Chemotherapy',
        'E': 'None of the above'
    }

    # Analysis based on the provided text
    # The text describes the treatment sequence for HER2+ Metastatic Breast Cancer (MBC).
    # First-line therapy is specified as THP (trastuzumab, pertuzumab, and a taxane).
    # The question asks for the treatment upon progression, which is the second-line therapy.
    
    key_sentence_1 = "Upon disease progression, TDXd is recommended in the second-line setting."
    key_sentence_2 = "In 2022, TDXd replaced T-DM1 as second-line therapy on the basis of superior results reported in the DESTINY-Breast03 trial."
    
    # TDXd is the abbreviation for fam-trastuzumab deruxtecan.
    correct_option = 'C'

    print(f"Question: {question}\n")
    print("Analysis of the text reveals the recommended treatment sequence from the NCCN.")
    print("-" * 70)
    print("Step 1: Identify the recommended first-line therapy.")
    print("The text states that a combination including a taxane and trastuzumab (as part of THP) is the preferred first-line regimen.\n")
    
    print("Step 2: Identify the recommended second-line therapy after progression.")
    print("The text directly addresses this:\n")
    print(textwrap.fill(f"'{key_sentence_1}'", width=70))
    print(textwrap.fill(f"'{key_sentence_2}'", width=70))
    print("\nTDXd stands for Fam-trastuzumab deruxtecan.")
    
    # Fulfilling the "output each number" requirement by highlighting key numerical/named entities
    trial_name = "DESTINY-Breast03"
    treatment_line = 2
    
    print("-" * 70)
    print("Final Answer Derivation:")
    print(f"Treatment Line: {treatment_line}nd line")
    print(f"Recommended Drug: {options[correct_option]} (TDXd)")
    print(f"Supporting Trial: {trial_name}")
    print(f"\nTherefore, the correct answer is C.")

solve_medical_question()