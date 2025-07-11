import textwrap

def solve_medical_question():
    """
    This script analyzes the provided text to answer the multiple-choice question.
    """

    # The question and answer choices from the user prompt.
    question = "Which of these has a category 1 recommendation from the NCCN for use in patients with HER2+ MBC whose disease progresses after treatment with trastuzumab and a taxane?"
    options = {
        "A": "Tucatinib trastuzumab and capecitabine",
        "B": "Trastuzumab emtansine",
        "C": "Fam-trastuzumab deruxtecan",
        "D": "Chemotherapy",
        "E": "None of the above"
    }

    # The correct answer key determined by analyzing the text.
    correct_answer_key = "C"

    # Explanation based on the provided text.
    explanation_steps = [
        "1. The question asks for the recommended therapy after the disease progresses on a first-line treatment that includes a taxane and trastuzumab for HER2+ MBC.",
        "2. The provided text states that 'Upon disease progression TDXd is recommended in the second-line setting.'",
        "3. The text also clarifies that 'In 2022 TDXd replaced T-DM1 as second-line therapy on the basis of superior results...'. T-DM1 corresponds to option B, trastuzumab emtansine.",
        "4. TDXd is the abbreviation for 'fam-trastuzumab deruxtecan', which corresponds to option C.",
        "5. For option A, the text states 'NCCN guidelines recommend tucatinib, trastuzumab, and capecitabine as the preferred regimen' for third-line therapy, not second-line."
    ]

    print("Analyzing the provided text to find the NCCN-recommended second-line therapy for HER2+ MBC:")
    print("-" * 70)
    for step in explanation_steps:
        # Wrap text for better readability in the terminal
        print("\n".join(textwrap.wrap(step, width=70)))

    print("-" * 70)
    print(f"Conclusion: The correct option is C, {options[correct_answer_key]}.")

    # Final answer in the required format.
    print(f"\n<<<{correct_answer_key}>>>")

# Execute the function to find and print the answer.
solve_medical_question()