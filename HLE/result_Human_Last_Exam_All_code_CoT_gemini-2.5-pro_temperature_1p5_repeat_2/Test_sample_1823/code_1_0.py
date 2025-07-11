def find_diagnosis():
    """
    This script analyzes the clinical information provided to determine the most likely diagnosis.
    """
    # Clinical Data
    patient_age = 1  # in years
    physical_findings = ["hypertrophic scarring", "erythema", "spasticity"]
    lab_results = "Negative anti-Mi-2"

    print("Patient Presentation Analysis:")
    print(f"- Age: {patient_age}-year-old")
    print(f"- Key Findings: {', '.join(physical_findings)}")
    print(f"- Lab Results: {lab_results}")
    print("-" * 30)
    print("Evaluating the Most Likely Diagnosis:")
    print(
        "The patient presents with symptoms involving both the skin (erythema, scarring) "
        "and muscles (spasticity)."
    )
    print(
        "We need a diagnosis that explains this combination of findings.\n"
    )
    print("Dermatomyositis is an inflammatory condition affecting the skin ('dermato-') and muscles ('-myositis').")
    print("- Erythema is a classic skin manifestation.")
    print("- Muscle inflammation causes weakness, which can lead to stiffness and contractures perceived as spasticity.")
    print("- Skin ulcerations from the disease can heal with hypertrophic scarring.")
    print("- The anti-Mi-2 antibody test, while specific, is often negative, especially in juvenile cases. "
          "A negative result does not exclude the diagnosis.")
    print("\nOther options like Ectropion, McArdle disease, McCune Albright syndrome, and Cataracts do not fit the "
          "full clinical picture of simultaneous skin and muscle pathology.")
    print("-" * 30)
    print("Conclusion: The combination of findings strongly points to Dermatomyositis.")

find_diagnosis()
<<<C>>>