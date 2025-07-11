def find_recommended_therapy():
    """
    This function determines the recommended therapy based on the provided text.
    The user is asking for the NCCN Category 1 recommended therapy for HER2+ MBC
    after the disease has progressed on a first-line regimen (trastuzumab and a taxane).
    This means we are looking for the second-line therapy.
    """

    # Information extracted from the provided text
    nccn_guidelines = {
        "first_line": "Trastuzumab, Pertuzumab, and a taxane (THP)",
        "second_line": "Fam-trastuzumab deruxtecan (TDXd)",
        "third_line": "Tucatinib, trastuzumab, and capecitabine"
    }

    answer_choices = {
        'A': 'Tucatinib trastuzumab and capecitabine',
        'B': 'Trastuzumab emtansine',
        'C': 'Fam-trastuzumab deruxtecan',
        'D': 'Chemotherapy',
        'E': 'None of the above'
    }

    # The question asks for the therapy after the first line, which is the second-line therapy.
    line_of_therapy_in_question = "second_line"
    recommended_drug = nccn_guidelines[line_of_therapy_in_question]

    # Find the matching answer choice
    correct_choice = None
    for choice, description in answer_choices.items():
        if description in recommended_drug:
            correct_choice = choice
            break

    print("The question asks for the NCCN-recommended therapy for HER2+ MBC after progression on first-line treatment.")
    print("According to the text, the recommended second-line therapy is Fam-trastuzumab deruxtecan (TDXd).")
    
    # Outputting the components of the answer
    print("\n--- Answer Components ---")
    print(f"Line of Therapy: 2")
    print(f"Recommended Drug: {recommended_drug}")
    print(f"Matching Answer Choice: {correct_choice}")
    print("-------------------------\n")


find_recommended_therapy()
print("<<<C>>>")