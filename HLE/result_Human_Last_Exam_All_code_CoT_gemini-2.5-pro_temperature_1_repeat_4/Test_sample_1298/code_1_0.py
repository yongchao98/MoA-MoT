import textwrap

def get_clinical_decision():
    """
    Analyzes the clinical case and determines the best course of action.
    """
    patient_profile = {
        "Diagnosis": "Metastatic Malignant Melanoma",
        "Complication": "Malignant pericardial effusion causing cardiac tamponade",
        "Key Finding": "Malignant cells found in pericardial fluid after pericardiocentesis."
    }

    reasoning = [
        "The patient has a confirmed systemic malignancy (metastatic melanoma).",
        "The immediate life-threatening issue (cardiac tamponade) has been stabilized by pericardiocentesis.",
        "The next priority is to treat the underlying widespread cancer.",
        "Symptomatic treatments (analgesia, diuretics) are insufficient as they do not address the root cause.",
        "Localized therapy (radiotherapy) is not the primary choice for systemic disease.",
        "Systemic therapy is required to target cancer cells throughout the body.",
        "Among the choices, Chemotherapy is the appropriate systemic treatment for the malignant cells."
    ]

    conclusion = "The next best step is to treat the underlying systemic cancer. Chemotherapy is a systemic therapy designed to kill malignant cells and is the most appropriate choice provided."
    
    best_choice = "D"

    print("Clinical Case Analysis:")
    print("-" * 30)
    print(f"Primary Diagnosis: {patient_profile['Diagnosis']}")
    print(f"Acute Complication: {patient_profile['Complication']}")
    print(f"Confirmatory Finding: {patient_profile['Key Finding']}")
    print("\nReasoning for Next Best Step:")
    for point in reasoning:
        print(f"- {point}")
    print("\n" + textwrap.fill(conclusion, width=70))
    print(f"\nTherefore, the best answer choice is: {best_choice}")


if __name__ == "__main__":
    get_clinical_decision()