def find_diagnosis():
    """
    Analyzes the patient's symptoms to determine the most likely diagnosis from the given options.
    """
    patient_age = 1  # year old
    symptoms = ["hypertrophic scarring", "erythema", "spasticity"]
    labs = "negative for anti-Mi-2"

    print("Patient Profile:")
    print(f"- Age: {patient_age} year old")
    print(f"- Symptoms: {', '.join(symptoms)}")
    print(f"- Labs: {labs}\n")

    print("Analysis of Answer Choices:")
    
    print("\nA. Ectropion: This is a condition where the eyelid turns outward. It does not explain the systemic symptoms of erythema and spasticity.")
    
    print("\nB. McArdle disease: This is a metabolic disorder affecting muscles. It typically causes exercise intolerance and muscle cramps, not the prominent skin findings (erythema, scarring) seen in this patient.")
    
    print("\nC. Dermatomyositis: This is an inflammatory disease affecting the skin and muscles.")
    print("  - Erythema (skin redness) is a hallmark sign.")
    print("  - Muscle involvement is key. While weakness is more common, inflammation and contractures can cause stiffness that may be perceived as spasticity.")
    print("  - In Juvenile Dermatomyositis, a common complication is calcinosis cutis (calcium deposits), which can ulcerate and lead to severe scarring, explaining the hypertrophic scars.")
    print("  - The anti-Mi-2 antibody is specific but not sensitive; a negative result does not rule out dermatomyositis.")
    print("  - Conclusion: This diagnosis aligns well with the combination of skin and muscle symptoms.")

    print("\nD. McCune Albright syndrome: This genetic disorder has a classic triad of fibrous dysplasia of bone, café-au-lait spots, and endocrine problems. This does not match the patient's symptoms.")

    print("\nE. Cataracts: This is a clouding of the eye's lens and does not account for the skin or muscle findings.")

    print("\nFinal Conclusion:")
    print("Based on the analysis, Dermatomyositis is the only option that can account for the combination of characteristic skin inflammation (erythema), significant muscle involvement (spasticity), and potential for severe scarring. It is the most likely diagnosis among the choices.")

if __name__ == "__main__":
    find_diagnosis()