def solve_diagnosis():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    # 1. Define the patient's key clinical features.
    patient_age = "1-year-old"
    physical_findings = {"Hypertrophic scarring", "Erythema", "Spasticity"}
    lab_results = "Negative for anti-Mi-2"

    print("Analyzing the clinical case based on the following key features:")
    print(f"- Patient Age: {patient_age}")
    print(f"- Physical Exam: {', '.join(sorted(list(physical_findings)))}")
    print(f"- Lab Results: {lab_results}\n")

    # 2. Evaluate each possible diagnosis.
    print("--- Evaluation of Answer Choices ---")

    # Choice A: Ectropion
    print("\nA. Ectropion:")
    print("   - This is an outward turning of the eyelid. It is a localized ocular finding and does not explain systemic issues like erythema on the body, hypertrophic scarring, or spasticity.")
    print("   - Ruling: Inconsistent with the patient's multi-system presentation.")

    # Choice B: McArdle disease
    print("\nB. McArdle disease:")
    print("   - This is a muscle glycogen storage disease causing exercise intolerance and muscle cramps, typically with onset in adolescence.")
    print("   - It does not cause the skin findings of erythema or hypertrophic scarring.")
    print("   - Ruling: Poor fit for the age and skin manifestations.")

    # Choice C: Dermatomyositis
    print("\nC. Dermatomyositis:")
    print("   - This is an inflammatory disease affecting both the skin and muscles, which aligns with the patient's symptoms.")
    print("   - Erythema is a classic skin finding (e.g., heliotrope rash).")
    print("   - In Juvenile Dermatomyositis (JDM), complications like calcinosis cutis can lead to skin ulceration and subsequent hypertrophic scarring.")
    print("   - Muscle involvement is key. While weakness is more typical than spasticity, 'spasticity' points towards a significant muscle pathology.")
    print(f"   - Crucially, the '{lab_results}' result is common in JDM (present in <10% of cases), making this diagnosis *more* likely in a {patient_age} patient, not less.")
    print("   - Ruling: Strongest candidate as it connects all key features.")

    # Choice D: McCune Albright syndrome
    print("\nD. McCune Albright syndrome:")
    print("   - This syndrome has a classic triad of fibrous dysplasia (bone disease), café-au-lait spots (a specific skin pigmentation), and endocrine problems. These do not match the patient's findings.")
    print("   - Ruling: Inconsistent with the presenting signs.")

    # Choice E: Cataracts
    print("\nE. Cataracts:")
    print("   - This is a clouding of the eye's lens. Like ectropion, it does not explain the skin and muscular symptoms.")
    print("   - Ruling: Inconsistent with the patient's multi-system presentation.")

    # 3. Formulate the final conclusion as a logical equation.
    print("\n--- Final Conclusion ---")
    print("The final diagnosis is determined by combining all the evidence:")
    final_equation = f"Patient Age ({patient_age}) + Skin Findings (Erythema, Hypertrophic scarring) + Muscle Findings (Spasticity) + Lab Result ({lab_results}) = Most Likely Diagnosis: Dermatomyositis (Juvenile Form)"
    print(final_equation)


solve_diagnosis()