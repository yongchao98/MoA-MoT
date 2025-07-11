def diagnose_patient():
    """
    Analyzes clinical findings to determine the most likely diagnosis from a set of choices.
    """
    # Patient Data
    patient_age = 1  # in years
    symptoms = ["hypertrophic scarring", "erythema", "spasticity"]
    labs = {"anti-Mi-2": "negative"}

    # Answer Choices
    choices = {
        "A": "Ectropion",
        "B": "McArdle disease",
        "C": "Dermatomyositis",
        "D": "McCune Albright syndrome",
        "E": "Cataracts"
    }

    print("Patient Profile:")
    print(f"- Age: {patient_age} year(s) old")
    print(f"- Symptoms: {', '.join(symptoms)}")
    print(f"- Key Lab Finding: anti-Mi-2 is {labs['anti-Mi-2']}")
    print("\n--- Diagnostic Reasoning ---")

    # Reasoning for excluding other options
    print("Excluding A & E: Ectropion and Cataracts are localized eye conditions and do not explain the systemic symptoms.")
    print("Excluding B: McArdle disease (infantile form) presents with hypotonia (low muscle tone), the opposite of spasticity.")
    print("Excluding D: McCune Albright syndrome presents with a different set of symptoms (bone dysplasia, specific skin spots, endocrine issues).")
    
    # Reasoning for the most likely diagnosis
    print("\nEvaluating Choice C: Dermatomyositis")
    print("- The patient's age (1) is consistent with Juvenile Dermatomyositis (JDM).")
    print("- Erythema (skin inflammation) is a classic feature.")
    print("- Hypertrophic scarring can result from severe skin vasculitis and ulceration seen in JDM.")
    print("- A negative anti-Mi-2 test is common in JDM.")
    print("- Spasticity, while not typical, can be a rare complication of CNS vasculitis in severe JDM.")
    
    most_likely_diagnosis_code = "C"
    most_likely_diagnosis_name = choices[most_likely_diagnosis_code]

    print("\n--- Conclusion ---")
    print(f"The patient's combination of symptoms and lab results makes '{most_likely_diagnosis_name}' the most likely diagnosis.")
    print(f"Final Answer Code: {most_likely_diagnosis_code}")

if __name__ == "__main__":
    diagnose_patient()