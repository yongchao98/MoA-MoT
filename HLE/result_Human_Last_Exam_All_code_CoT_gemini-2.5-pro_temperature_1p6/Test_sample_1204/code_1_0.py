import textwrap

def explain_treatment_plan():
    """
    This function outlines the reasoning for the selected treatment options.
    """
    analysis = {
        'I. Counsel patient on stopping cannabis.': 'This is a critical first step. The patient\'s heavy cannabis use is a significant confounding variable that likely worsens his anxiety and is a primary cause of his severe insomnia. It\'s impossible to accurately assess his underlying conditions or the effectiveness of his medications while this continues.',
        'III. Order a urine drug test.': 'This is essential for patient safety and effective treatment planning. Given his history of polysubstance use, a UDS provides objective data to confirm his cannabis use and screen for the return of other substances like cocaine or the use of unmentioned drugs.',
        'IV. Prescribe melatonin for insomnia.': 'The patient is struggling mightily with insomnia. While his cannabis use is the most likely cause, providing a safe, non-addictive sleep aid like melatonin directly addresses a chief complaint. This can provide some relief and help build the therapeutic alliance needed to tackle the more difficult issue of cannabis cessation.',
        'Why other options are less prioritized:': {
            'II. Hospitalization': 'This is a drastic measure and not an appropriate initial step in an outpatient setting. It should be considered only if outpatient management fails.',
            'V. Change AUD meds': 'His AUD is in remission. It is not wise to alter a medication regimen that is currently associated with successful remission without a compelling reason.',
            'VI. Start atomoxetine': 'The patient is already taking this medication. Furthermore, addressing questionable ADHD is a lower priority than managing his substance use, insomnia, and mood instability.'
        }
    }

    print("The three treatment options that should be prioritized immediately are I, III, and IV.\n")
    print("Here is the reasoning for each choice:\n")

    for option, reason in analysis.items():
        if isinstance(reason, str):
            print(f"Chosen Option: {option}")
            print(textwrap.fill(f"Reason: {reason}", width=80))
            print("-" * 80)
        else:
            print(f"{option}\n")
            for sub_option, sub_reason in reason.items():
                print(f"   - {sub_option}")
                print(textwrap.fill(f"     {sub_reason}", width=80))
            print("-" * 80)

explain_treatment_plan()

# The final answer is the letter corresponding to the combination I, III, and IV.
# A=I,II,III; B=I,II,IV; C=I,II,V; D=I,II,VI; E=II,III,IV; F=II,III,V;
# G=II,III,VI; H=II,IV,V; I=II,IV,VI; J=III,IV,V; K=III,IV,VI; L=I,III,IV
print("\nFinal Answer Choice Code:")
print("<<<L>>>")