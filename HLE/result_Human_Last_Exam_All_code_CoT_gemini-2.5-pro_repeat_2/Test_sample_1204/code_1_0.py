def solve_clinical_case():
    """
    Analyzes the clinical case and determines the three most appropriate immediate treatment options.
    """
    options = {
        'I': "Counsel patient on stopping cannabis.",
        'II': "Ask patient to request admission to the hospital so he can be detoxed off all of his psych meds.",
        'III': "Order a urine drug test.",
        'IV': "Prescribe melatonin for insomnia.",
        'V': "Discontinue acamprosate and increase dosage of naltrexone for AUD.",
        'VI': "Start atomoxetine."
    }

    # The chosen options are I, III, and IV.
    chosen_keys = ['I', 'III', 'IV']
    
    print("The three treatment options that should be prioritized immediately are:\n")

    # Print and justify the first choice
    key1 = chosen_keys[0]
    print(f"{key1}. {options[key1]}")
    print("   - Rationale: The patient's heavy daily cannabis use is a major confounding factor for his anxiety, mood, and insomnia. It is a substance use disorder in itself that must be addressed to allow for accurate assessment and effective treatment of his other conditions.\n")

    # Print and justify the second choice
    key2 = chosen_keys[1]
    print(f"{key2}. {options[key2]}")
    print("   - Rationale: A urine drug test provides essential objective data on the patient's substance use. Given his history and his request for stimulants (which have abuse potential), it is critical to confirm reported use and screen for other unreported substances to ensure safety and guide treatment.\n")

    # Print and justify the third choice
    key3 = chosen_keys[2]
    print(f"{key3}. {options[key3]}")
    print("   - Rationale: Insomnia is a major complaint causing the patient significant distress. While likely multifactorial, prescribing a safe, non-addictive agent like melatonin addresses a chief complaint directly. This can provide symptomatic relief and help build the therapeutic alliance needed to tackle the more complex issues like cannabis use and polypharmacy.\n")

solve_clinical_case()
print("<<<L>>>")