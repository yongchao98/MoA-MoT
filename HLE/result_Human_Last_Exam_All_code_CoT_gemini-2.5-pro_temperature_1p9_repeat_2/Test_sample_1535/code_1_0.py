def solve_clinical_case():
    """
    This function analyzes the patient's symptoms and signs to determine the expected location of a rash.
    """
    # Key information from the case presentation
    patient_history = [
        "muscle weakness",
        "myalgia (muscle pain)",
        "arthralgia (joint pain)",
        "fatigue"
    ]
    key_physical_finding = "periorbital erythema (redness around the eyes)"

    print("Step 1: Analyze the patient's key symptoms and signs.")
    print("The patient exhibits systemic symptoms suggestive of an inflammatory myopathy:")
    for symptom in patient_history:
        print(f"- {symptom}")
    print("\nThe most specific physical finding is:")
    print(f"- {key_physical_finding}\n")

    print("Step 2: Formulate a likely diagnosis based on the evidence.")
    print("The combination of muscle inflammation (myositis) and skin findings (dermatitis) strongly suggests a diagnosis of Dermatomyositis.")
    print("The finding of 'periorbital erythema' is a classic description of a 'Heliotrope rash,' which is a hallmark sign of this condition.\n")

    print("Step 3: Relate the diagnosis to the anatomical location in question.")
    print("A Heliotrope rash is a violaceous erythematous rash located specifically on the upper eyelids, sometimes extending around the eyes.\n")

    print("Step 4: Evaluate the answer choices based on this conclusion.")
    print("A. Dorsum of the hands (Gottron's papules - also in dermatomyositis, but not the sign described)")
    print("B. Nose (Malar rash of Lupus)")
    print("C. Eyelids (The exact location of a Heliotrope rash, matching 'periorbital erythema')")
    print("D. Groin (Not a characteristic location)")
    print("E. Shoulders (Shawl sign - also in dermatomyositis, but not the sign described)\n")
    
    print("Conclusion: The physical finding of periorbital erythema directly corresponds to a rash on the eyelids.")

# Execute the analysis
solve_clinical_case()
<<<C>>>