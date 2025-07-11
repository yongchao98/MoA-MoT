import textwrap

def analyze_clinical_case():
    """
    Analyzes the provided clinical vignette to determine the best categorization for the patient's pathology.
    """
    case_summary = {
        "Age": "60 years old",
        "Chief Complaint": "Memory loss",
        "Daughter's Report": [
            "Forgets to feed himself",
            "Weight loss",
            "Disorientation to day, month, or year"
        ],
        "Patient's Statements": [
            "Correctly recalls 3 objects immediately",
            "Denies daughter's claims (denial/lack of insight)",
            "Claims weight loss is due to a 'rare tapeworm' (confabulation)"
        ],
        "Medical History": "Chronic venous insufficiency",
        "Pertinent Negatives": ["Hypertension", "Cirrhosis"],
        "Psychosocial History": "10 pack years of smoking",
        "Physical Exam": "Normal"
    }

    print("--- Analyzing the Clinical Case ---")
    print("Patient presents with significant memory loss, disorientation, and confabulation (fabricating an explanation for his weight loss). Let's evaluate the options:")
    print("-" * 35)

    # Analysis of each option
    analysis = {
        "A": "Short-term memory: This aligns directly with the core symptoms. The patient forgets to eat and is disoriented to time, both classic signs of severe short-term memory impairment. His statement about a 'tapeworm' is a confabulation, a key feature in some amnestic syndromes where patients invent stories to fill memory gaps. This is the most fitting category.",
        "B": "Restrictive cardiomyopathy: This is a heart condition. There are no signs or symptoms mentioned (like shortness of breath, chest pain, or abnormal heart sounds) and the physical exam is normal. This is incorrect.",
        "C": "Hepatic encephalopathy: This is confusion caused by severe liver failure. The case explicitly states the patient does *not* have cirrhosis, making this diagnosis highly unlikely.",
        "D": "Parasitic infection: The patient's claim of a 'rare tapeworm' is presented in a context that strongly suggests it is a confabulation to explain his weight loss, which is actually due to him forgetting to eat. There is no objective evidence of an infection.",
        "E": "ATP depletion: This is a cellular-level process, not a clinical diagnosis or a category of pathology. While brain cell dysfunction involves energy issues, it is too general and not how a clinician would categorize this patient's condition."
    }

    for option, reason in analysis.items():
        # Use textwrap to format the output nicely
        wrapped_text = textwrap.fill(f"Option {option}: {reason}", width=80)
        print(wrapped_text)
        print("-" * 35)

    print("\nConclusion: The patient's entire clinical picture, from memory loss to disorientation and confabulation, points to a primary deficit in short-term memory.")


if __name__ == "__main__":
    analyze_clinical_case()