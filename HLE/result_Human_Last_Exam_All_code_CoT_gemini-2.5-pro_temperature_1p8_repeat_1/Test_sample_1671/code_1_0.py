import textwrap

def find_mbc_treatment_recommendation():
    """
    This script analyzes the treatment guidelines for HER2+ Metastatic Breast Cancer (MBC)
    as described in the provided text to answer the user's question.
    """

    # Information extracted from the provided text about NCCN guidelines
    nccn_recommendations = {
        "first-line": "THP (chemotherapy with trastuzumab and pertuzumab)",
        "second-line": "Fam-trastuzumab deruxtecan (TDXd)",
        "third-line": "Tucatinib, trastuzumab, and capecitabine",
        "historical_second_line": "Trastuzumab emtansine (T-DM1)"
    }

    # The user's question and the available choices
    question = "Which of these has a category 1 recommendation from the NCCN for use in patients with HER2+ MBC whose disease progresses after treatment with trastuzumab and a taxane?"
    
    choices = {
        "A": "Tucatinib trastuzumab and capecitabine",
        "B": "Trastuzumab emtansine",
        "C": "Fam-trastuzumab deruxtecan",
        "D": "Chemotherapy",
        "E": "None of the above"
    }

    print("Question:")
    print(textwrap.fill(question, width=80))
    print("\nAnswer Choices:")
    for key, value in choices.items():
        print(f"{key}. {value}")
    print("-" * 30)

    # The question asks for the therapy after first-line treatment, which is the second-line therapy.
    required_therapy_line = "second-line"
    correct_treatment = nccn_recommendations[required_therapy_line]
    
    # Determine the correct answer choice
    correct_choice_letter = ""
    for letter, choice_text in choices.items():
        # Check if the recommended treatment name is contained within the choice text
        if choice_text.lower() in correct_treatment.lower():
            correct_choice_letter = letter
            break

    # Print the reasoning based on the text
    print("\nReasoning:")
    print(f"The question asks for the standard of care after first-line therapy, which is the '{required_therapy_line}' treatment.")
    print(f"According to the provided text, the NCCN recommendation for the {required_therapy_line} setting is '{correct_treatment}'.")
    print(f"The text specifically states: 'Upon disease progression TDXd is recommended in the second-line setting.'")
    print(f"It also clarifies: '...fam-trastuzumab deruxtecan (TDXd) are approved for HER2+ MBC.'")
    print(f"This matches option {correct_choice_letter}.")
    
    # Final answer in the required format
    print(f"\nFinal Answer:")
    print(f"<<<{correct_choice_letter}>>>")

# Execute the function
find_mbc_treatment_recommendation()