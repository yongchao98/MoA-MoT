import textwrap

def solve_clinical_case():
    """
    Analyzes the clinical vignette and determines the best course of action.
    """

    # Step 1 & 2: Analyze presentation and formulate diagnosis
    patient_findings = {
        "History": "20-year-old woman with persistent fatigue, weight loss.",
        "Skin": "Changing mole (darkening, growing, irregular border, raised) and new dark spots on arms/chest.",
        "Systemic Symptoms": "Dull right hip ache, shortness of breath, chest discomfort, swelling in legs/abdomen.",
        "Physical Exam": "Muffled heart sounds, jugular venous distention (JVD).",
        "Investigations": "Pericardiocentesis fluid shows malignant cells, high protein, and high LDH."
    }

    diagnosis_reasoning = """
    The patient's history of a changing mole is highly suspicious for malignant melanoma.
    The new dark spots, hip ache, and systemic symptoms suggest metastatic disease.
    The muffled heart sounds, JVD, and malignant cells in the pericardial fluid confirm a malignant pericardial effusion causing cardiac tamponade, a life-threatening complication of metastatic cancer.
    The primary diagnosis is Metastatic Malignant Melanoma.
    """

    # Step 3: Identify Acute vs. Chronic Issues
    management_context = """
    The immediate life-threatening condition, cardiac tamponade, has been managed with pericardiocentesis.
    The next best step must address the underlying cause of all these problems: the widespread metastatic cancer.
    """

    # Step 4: Evaluate Management Options
    options = {
        "A": "Prescribe meloxicam to manage the persistent fatigue. (Symptomatic relief, does not treat cancer)",
        "B": "Prescribe low-dose analgesia. (Symptomatic relief, does not treat cancer)",
        "C": "None of the choices.",
        "D": "Chemotherapy to kill the malignant cells. (Systemic treatment that targets cancer cells throughout the body)",
        "E": "Immunosuppression. (Contraindicated; immunotherapy, which boosts the immune system, is used for melanoma)",
        "F": "Rapid antibiotic infusion. (No evidence of infection)",
        "G": "Radiotherapy to treat the malignant cells. (Local treatment, not suitable for widespread disease)",
        "H": "Diuretics to reduce the fluid overload. (Symptomatic relief, does not treat the cause of fluid accumulation)"
    }

    # Step 5: Select the Best Answer
    conclusion = """
    The patient has a systemic malignancy (metastatic melanoma). Therefore, a systemic treatment is required.
    Of the choices provided, Chemotherapy is a systemic therapy designed to kill malignant cells throughout the body.
    While modern melanoma treatment often involves immunotherapy or targeted therapy, chemotherapy is the correct category of treatment listed to address the widespread disease. It is the most appropriate next step to control the cancer and prevent further complications like the re-accumulation of pericardial fluid.
    """

    correct_choice = "D"

    print("Clinical Case Analysis:")
    print("=" * 25)
    print("\n[Patient Presentation]")
    for key, value in patient_findings.items():
        print(f"- {key}: {value}")

    print("\n[Diagnosis]")
    print(textwrap.fill(diagnosis_reasoning.strip(), width=80))

    print("\n[Management Context]")
    print(textwrap.fill(management_context.strip(), width=80))

    print("\n[Evaluation of Choices]")
    for key, value in options.items():
        print(f"- {key}: {value}")

    print("\n[Conclusion]")
    print(textwrap.fill(conclusion.strip(), width=80))

    print("\n" + "="*25)
    print(f"The final answer is option {correct_choice}.")


solve_clinical_case()

# The final answer is derived from the clinical reasoning above.
# The patient has metastatic melanoma. The immediate complication (cardiac tamponade) has been stabilized.
# The next step is systemic treatment for the cancer. Chemotherapy is the only systemic anti-cancer therapy listed.
print("<<<D>>>")