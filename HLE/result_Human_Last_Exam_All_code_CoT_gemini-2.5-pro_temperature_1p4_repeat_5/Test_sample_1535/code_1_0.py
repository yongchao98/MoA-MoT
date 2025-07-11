def diagnose_rash_location():
    """
    Analyzes a patient's clinical data to determine the expected location of a rash.
    """
    # Key information from the case
    patient_age = 45
    key_finding = "periorbital erythema"
    supporting_findings = ["myalgia", "arthralgia", "muscle weakness"]
    
    # Define the reasoning process as a series of steps
    print("Clinical Reasoning Process:")
    print("Step 1: The most specific physical exam finding is '{}'.".format(key_finding))
    print("Step 2: This finding is characteristic of a sign known as a 'Heliotrope rash'.")
    print("Step 3: A Heliotrope rash, especially in a {}-year-old patient with {}, is a classic sign of Dermatomyositis.".format(patient_age, ", ".join(supporting_findings)))
    print("Step 4: The specific anatomical location of a Heliotrope rash is the eyelids.")
    
    # Summarize the logical flow in an equation format, including numbers as requested
    print("\n---")
    print("Conclusion represented as a logical equation:")
    print("The final answer is derived from the following steps:")
    # The following line prints the numbers 1, 2, 3, and 4 in the final equation.
    print("Step 1 (Symptom) + Step 2 (Sign) + Step 3 (Diagnosis) leads to Step 4 (Location: Eyelids)")

# Execute the diagnostic function
diagnose_rash_location()