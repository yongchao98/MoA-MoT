import textwrap

def explain_diagnosis():
    """
    Explains the reasoning for choosing the best next diagnostic step based on the clinical case.
    """
    case_summary = {
        "Diagnosis": "The clinical presentation is highly suggestive of allergic contact dermatitis, likely from new workout clothes (textile dermatitis).",
        "Key Finding": "The rash spares the axillary vault and is located on the periphery where clothing causes friction, pointing away from deodorant and towards clothing as the cause."
    }

    options = {
        'A': "Skin biopsy is invasive and not the first-line test to identify a specific allergen. It's used when the diagnosis is unclear.",
        'B': "KOH preparation is used to diagnose fungal infections. The presentation is more typical of contact dermatitis, not a fungal infection.",
        'C': "Topical steroid is a treatment, not a diagnostic test. While its effectiveness can support an inflammatory diagnosis, it doesn't identify the cause.",
        'D': "Patch test is the gold standard for diagnosing allergic contact dermatitis. It can identify the specific chemical (allergen) in the clothing that is causing the reaction. This is the most logical next step to confirm the suspected diagnosis.",
        'E': "None of the above is incorrect because the patch test is a very appropriate next step."
    }

    print("Analyzing the Clinical Case and Diagnostic Options:")
    print("-" * 50)
    for key, value in case_summary.items():
        print(f"{key}: {value}")
    
    print("\nEvaluating the Answer Choices:")
    print("-" * 50)
    for option, explanation in options.items():
        # Wrap text for better readability in terminal
        wrapped_explanation = textwrap.fill(f"{option}. {explanation}", width=80)
        print(wrapped_explanation)
    
    conclusion = "\nConclusion: The best next step to confirm the diagnosis of allergic contact dermatitis and identify the specific trigger is the patch test."
    print(conclusion)

explain_diagnosis()

print("\n<<<D>>>")