import sys

def solve_clinical_vignette():
    """
    Analyzes the clinical case and determines the three most prioritized immediate treatment options.
    """
    
    # Step 1: Analyze the core problems presented in the vignette.
    analysis = {
        "Problem 1: Substance Use": "Patient has an active Cannabis Use Disorder, which is likely worsening his insomnia and anxiety and reducing the effectiveness of his medications.",
        "Problem 2: Polypharmacy": "The patient is on a complex and potentially problematic combination of medications (e.g., two antidepressants).",
        "Problem 3: Diagnostic Clarity": "There is a history of multiple psychiatric and substance use disorders, and the patient is requesting stimulants for 'questionable' ADHD.",
        "Primary Complaint": "The patient's chief complaints are severe insomnia and anxiety."
    }

    # Step 2: Evaluate each potential action.
    evaluation = {
        "I. Counsel patient on stopping cannabis.": "HIGHEST PRIORITY. This addresses the most likely root cause of his current primary complaints (insomnia, anxiety) and medication ineffectiveness.",
        "II. Ask patient to request admission to the hospital so he can be detoxed off all of his psych meds.": "INAPPROPRIATE. This is a drastic, high-risk action. Medication changes should be gradual and outpatient.",
        "III. Order a urine drug test.": "HIGHEST PRIORITY. This is essential for objective data gathering, confirming substance use, and ensuring patient safety before considering any medication changes, especially given his request for stimulants.",
        "IV. Prescribe melatonin for insomnia.": "REASONABLE PRIORITY. This is a safe, low-risk, non-addictive intervention that directly addresses a major patient complaint (insomnia) and can help build therapeutic rapport while the core issue of cannabis use is addressed.",
        "V. Discontinue acamprosate and increase dosage of naltrexone for AUD.": "LOW PRIORITY. There is no clinical indication to change treatment for a disorder that is in remission.",
        "VI. Start atomoxetine.": "INVALID. The patient is already taking atomoxetine according to the case description."
    }

    # Step 3: Select the top three priorities.
    # Based on the evaluation, the most logical and immediate priorities are I, III, and IV.
    prioritized_actions = ["I", "III", "IV"]
    
    print("Step-by-Step Reasoning:")
    print("1. The most immediate issue to address is the patient's heavy daily cannabis use, as it is a major contributing factor to his insomnia and anxiety. Therefore, counseling on cessation (I) is a top priority.")
    print("2. Given the patient's substance use history and his request for stimulants, obtaining objective data is crucial. A urine drug test (III) confirms active use and screens for other substances, which is essential for safe and effective treatment planning.")
    print("3. While addressing the root cause (cannabis), it is also important to address the patient's severe distress from insomnia. Prescribing melatonin (IV) is a safe, non-addictive, low-risk option to provide symptomatic relief.")
    print("\nConclusion: The three most appropriate immediate actions are to counsel the patient on stopping cannabis, order a urine drug test, and prescribe melatonin for insomnia.")
    
    final_answer_choice = "L"
    
    # The prompt asks to output the final equation, which in this context means the selected options.
    print(f"\nThe selected combination of options is: {prioritized_actions[0]}, {prioritized_actions[1]}, {prioritized_actions[2]}")
    
    # Final answer in the required format
    sys.stdout.write(f'<<<{final_answer_choice}>>>')

solve_clinical_vignette()