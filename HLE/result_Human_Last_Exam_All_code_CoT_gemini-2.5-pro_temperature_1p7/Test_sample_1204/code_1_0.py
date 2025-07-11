import collections

def solve_clinical_case():
    """
    Analyzes a clinical vignette to determine the three most appropriate
    immediate treatment options from a given list.
    """

    # Define the treatment options
    options = {
        "I": "Counsel patient on stopping cannabis.",
        "II": "Ask patient to request admission to the hospital so he can be detoxed off all of his psych meds.",
        "III": "Order a urine drug test.",
        "IV": "Prescribe melatonin for insomnia.",
        "V": "Discontinue acamprosate and increase dosage of naltrexone for AUD.",
        "VI": "Start atomoxetine."
    }

    # Step 1: Analyze each option's priority and provide a rationale.
    # We will score them: 3 = High Priority, 2 = Medium Priority, 1 = Low Priority, 0 = Inappropriate/Invalid.
    analysis = collections.OrderedDict()

    analysis["I"] = {
        "text": options["I"],
        "priority": 3,
        "rationale": "High Priority. The patient's heavy daily cannabis use is a major contributing factor to his insomnia and poorly controlled anxiety. Addressing the active substance use disorder is a critical first step."
    }
    analysis["III"] = {
        "text": options["III"],
        "priority": 3,
        "rationale": "High Priority. A urine drug screen is essential to get objective data. It confirms the patient's reported cannabis use and screens for other potential substances, which is vital given his history of cocaine use disorder."
    }
    analysis["IV"] = {
        "text": options["IV"],
        "priority": 2,
        "rationale": "Medium-High Priority. The patient is struggling 'mightily with insomnia'. Melatonin is a safe, non-addictive, low-risk intervention that can provide immediate symptomatic relief and help build a therapeutic alliance while addressing the root causes."
    }
    analysis["V"] = {
        "text": options["V"],
        "priority": 1,
        "rationale": "Low Priority. The patient's alcohol use disorder is in remission. There is no immediate clinical indication to change his maintenance medications (acamprosate, naltrexone). The focus should be on the acute, presenting problems."
    }
    analysis["II"] = {
        "text": options["II"],
        "priority": 0,
        "rationale": "Inappropriate. Abruptly stopping all psychotropic medications ('detox') is dangerous and could lead to severe withdrawal syndromes and significant clinical destabilization. Medication adjustments should be gradual and thoughtful."
    }
    analysis["VI"] = {
        "text": options["VI"],
        "priority": 0,
        "rationale": "Invalid. The clinical vignette states the patient is already taking atomoxetine 80 mg qD. Therefore, this action is redundant and incorrect."
    }

    print("Step-by-step analysis of treatment options:")
    for key, value in analysis.items():
        print(f"\nOption {key}: {value['text']}")
        print(f"Rationale: {value['rationale']}")

    # Step 2: Sort the options by priority to find the top three.
    sorted_options = sorted(analysis.items(), key=lambda item: item[1]['priority'], reverse=True)

    # Step 3: Select the top three options.
    top_three = sorted_options[:3]
    top_three_keys = sorted([item[0] for item in top_three]) # Sort keys alphabetically for consistent matching

    print("\n--------------------------------------------------")
    print("Conclusion of Analysis:")
    print("The three most critical immediate actions are:")
    for key, value in top_three:
        print(f"- Option {key}: {value['text']}")
    print("--------------------------------------------------")

    # Define all possible answer choices
    answer_choices = {
        "A": ["I", "II", "III"], "B": ["I", "II", "IV"], "C": ["I", "II", "V"], "D": ["I", "II", "VI"],
        "E": ["II", "III", "IV"], "F": ["II", "III", "V"], "G": ["II", "III", "VI"], "H": ["II", "IV", "V"],
        "I": ["II", "IV", "VI"], "J": ["III", "IV", "V"], "K": ["III", "IV", "VI"], "L": ["I", "III", "IV"],
        "M": ["I", "III", "V"], "N": ["I", "III", "VI"], "O": ["I", "IV", "V"], "P": ["I", "IV", "VI"],
        "Q": ["II", "V", "VI"], "R": ["III", "V", "VI"]
    }

    # Step 4: Find the matching answer choice.
    final_answer_letter = "Unknown"
    for letter, combo in answer_choices.items():
        if sorted(combo) == top_three_keys:
            final_answer_letter = letter
            break

    print(f"\nThe best combination of options is {', '.join(top_three_keys)}.")
    print(f"This corresponds to answer choice: {final_answer_letter}")

if __name__ == "__main__":
    solve_clinical_case()
