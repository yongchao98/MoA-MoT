import textwrap

def analyze_case():
    """
    Analyzes the clinical vignette to determine the best categorization of the patient's pathology.
    """
    vignette_summary = {
        "Chief Complaint": "Memory loss",
        "Key Symptoms": [
            "Forgets to feed himself, leading to weight loss",
            "Disoriented to day, month, and year",
            "Confabulation (fabricates a story about a 'rare tapeworm' to explain weight loss)",
            "Lack of insight (denies his daughter's accurate report)"
        ],
        "Pertinent History": "Chronic venous insufficiency, 10 pack-year smoking history",
        "Pertinent Negatives": "No hypertension, no cirrhosis",
        "Exam Findings": "Physical exam is normal. Can recall 3 objects immediately after being told."
    }

    print("Step 1: Analyzing the patient's presentation.")
    print("---------------------------------------------")
    for key, value in vignette_summary.items():
        if isinstance(value, list):
            print(f"- {key}:")
            for item in value:
                print(f"  - {item}")
        else:
            print(f"- {key}: {value}")
    print("\nThe core features are severe memory deficits, disorientation, and confabulation.\n")

    print("Step 2: Evaluating the answer choices.")
    print("---------------------------------------")

    # A. Short-term memory
    analysis_a = """
    A. Short-term memory: This directly describes the central problem. The patient's inability to recall recent events (like whether he has eaten or the current date) and his need to invent stories (confabulation) to fill in these memory gaps are hallmark signs of a severe pathology affecting the short-term memory system. This is the most accurate description of the primary functional deficit among the choices.
    """
    print(textwrap.dedent(analysis_a))

    # B. Restrictive cardiomyopathy
    analysis_b = """
    B. Restrictive cardiomyopathy: This is a heart condition. There are no signs or symptoms mentioned (e.g., shortness of breath, edema, abnormal heart sounds) to support this diagnosis. The physical exam is explicitly normal.
    """
    print(textwrap.dedent(analysis_b))

    # C. Hepatic encephalopathy
    analysis_c = """
    C. Hepatic encephalopathy: This is a brain dysfunction caused by severe liver failure. The vignette explicitly states 'no cirrhosis' as a pertinent negative, making this diagnosis highly unlikely.
    """
    print(textwrap.dedent(analysis_c))

    # D. Parasitic infection
    analysis_d = """
    D. Parasitic infection: The patient mentions a 'tapeworm', but this is a classic example of confabulation—a symptom of his memory disorder—not a real diagnosis. There is no objective evidence of an infection.
    """
    print(textwrap.dedent(analysis_d))

    # E. ATP depletion
    analysis_e = """
    E. ATP depletion: This is a non-specific biochemical process at the cellular level. It is not a clinical diagnosis or a useful way to categorize a patient's pathology in this context. It's too general and not clinically actionable.
    """
    print(textwrap.dedent(analysis_e))

    print("Step 3: Conclusion.")
    print("-------------------")
    conclusion = """
    The patient's clinical picture is dominated by a profound defect in forming and retrieving recent memories. The other options are either unsupported, contradicted by the provided information, or too non-specific. Therefore, categorizing the pathology as an issue of 'Short-term memory' is the most appropriate choice.
    """
    print(textwrap.dedent(conclusion))

if __name__ == "__main__":
    analyze_case()
    print("<<<A>>>")
