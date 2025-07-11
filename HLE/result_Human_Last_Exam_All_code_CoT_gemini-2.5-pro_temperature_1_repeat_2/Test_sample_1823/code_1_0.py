def solve_medical_case():
    """
    Analyzes the clinical vignette and determines the most likely diagnosis.
    """

    # Patient Data
    patient_age = 1  # years
    physical_findings = ["hypertrophic scarring", "erythema", "spasticity"]
    lab_results = "Negative for anti-Mi-2"

    # Print the analysis steps
    print("Step 1: Analyzing the patient's clinical presentation.")
    print(f"- A {patient_age}-year-old patient indicates a pediatric, likely juvenile, condition.")
    print(f"- The symptoms involve both skin ({physical_findings[0]}, {physical_findings[1]}) and muscle ({physical_findings[2]}).")
    print(f"- The lab result ({lab_results}) is a key piece of diagnostic information.")
    print("\nStep 2: Evaluating the differential diagnosis.")
    print("- A. Ectropion & E. Cataracts: Incorrect. These are eye conditions and do not explain systemic skin and muscle symptoms.")
    print("- B. McArdle disease: Incorrect. This typically presents later in life and does not cause these skin findings.")
    print("- D. McCune Albright syndrome: Incorrect. This syndrome has a different classic triad of symptoms (bone, skin spots, endocrine).")
    print("- C. Dermatomyositis: Correct. The juvenile form (JDM) fits the profile. It affects skin (erythema) and muscle (can cause contractures/spasticity). Skin ulceration from vasculopathy can lead to scarring. A negative anti-Mi-2 test is common in JDM.")

    print("\nStep 3: Final Conclusion.")
    final_answer = "C. Dermatomyositis"
    print(f"The combination of skin and muscle symptoms in a young child, especially with a negative anti-Mi-2 lab, strongly points to {final_answer}.")

# Execute the function to display the reasoning
solve_medical_case()