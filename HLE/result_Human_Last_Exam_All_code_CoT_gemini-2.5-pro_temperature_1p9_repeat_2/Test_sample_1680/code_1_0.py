import textwrap

def analyze_medical_case():
    """
    Analyzes the provided medical case to determine the best pathology categorization.
    """
    case_summary = {
        "Age": "60 years old",
        "Chief Complaint": "Memory loss",
        "Key Symptoms": [
            "Forgets to feed himself, leading to weight loss",
            "Disoriented to day, month, and year",
            "Confabulation (invents 'rare tapeworm' story to explain weight loss)",
            "Lacks insight into his condition"
        ],
        "Pertinent History": "10 pack years smoking",
        "Pertinent Negatives": "No hypertension, no cirrhosis"
    }

    answer_choices = {
        "A": "Short-term memory",
        "B": "Restrictive cardiomyopathy",
        "C": "Hepatic encephalopathy",
        "D": "Parasitic infection",
        "E": "ATP depletion"
    }

    print("Analyzing the medical case step-by-step:\n")

    print("1. Patient Presentation Analysis:")
    print(f"   - The patient (age {case_summary['Age']}) presents with classic signs of severe anterograde amnesia (inability to form new memories).")
    print(f"   - Evidence: He forgets to perform essential activities like eating and is disoriented to time.")
    print(f"   - A key feature is confabulation: He invents a plausible story ('rare tapeworm') to fill gaps in his memory, which points strongly to a primary memory disorder.\n")

    print("2. Evaluation of Answer Choices:")
    # Choice A
    reasoning_a = "This choice directly addresses the core clinical features. The patient's inability to recall recent events, disorientation, and confabulation are all hallmarks of a significant short-term memory pathology."
    print(f"   - Choice A ({answer_choices['A']}): CORRECT.")
    print(textwrap.fill(f"     Reasoning: {reasoning_a}\n", width=80, initial_indent='       ', subsequent_indent='       '))

    # Choice C
    reasoning_c = f"This is brain dysfunction due to liver failure. The case explicitly states a pertinent negative of 'no cirrhosis', making this diagnosis highly unlikely."
    print(f"   - Choice C ({answer_choices['C']}): INCORRECT.")
    print(textwrap.fill(f"     Reasoning: {reasoning_c}\n", width=80, initial_indent='       ', subsequent_indent='       '))

    # Choice D
    reasoning_d = "The 'parasitic infection' is not a diagnosis but a confabulation from the patient. It's a symptom of his memory disorder, not the cause of his condition. His weight loss is better explained by him forgetting to eat."
    print(f"   - Choice D ({answer_choices['D']}): INCORRECT.")
    print(textwrap.fill(f"     Reasoning: {reasoning_d}\n", width=80, initial_indent='       ', subsequent_indent='       '))

    # Other Choices
    reasoning_other = "Restrictive cardiomyopathy (B) is a heart condition with no supporting evidence. ATP depletion (E) is a cellular mechanism, not a clinical diagnosis for this syndrome."
    print("   - Choices B and E: INCORRECT.")
    print(textwrap.fill(f"     Reasoning: {reasoning_other}\n", width=80, initial_indent='       ', subsequent_indent='       '))

    print("Conclusion: The patient's entire clinical picture is best categorized as a pathology of short-term memory.")

# Run the analysis
analyze_medical_case()