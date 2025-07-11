def solve_medical_case():
    """
    This function analyzes the clinical vignette to identify the location of the expected rash.
    """

    # Key findings from the case description
    muscle_symptoms = ["muscle weakness", "myalgia", "arthralgia"]
    skin_finding = "periorbital erythema (redness around the eyes)"

    # Diagnosis based on the combination of key findings
    # The constellation of muscle weakness and characteristic skin findings strongly suggests Dermatomyositis.
    diagnosis = "Dermatomyositis"

    # Explanation of characteristic rashes in Dermatomyositis
    heliotrope_rash_description = "A violaceous (purplish-red) rash on the upper eyelids."
    
    print("Clinical Reasoning Steps:")
    print("1. The patient's symptoms of muscle weakness and myalgia are characteristic of an inflammatory myopathy.")
    print(f"2. The key dermatological finding is '{skin_finding}'.")
    print(f"3. This combination of muscle and skin symptoms is classic for a condition called {diagnosis}.")
    print(f"4. A hallmark sign of {diagnosis} is the 'Heliotrope rash,' which is described as: '{heliotrope_rash_description}'")
    print("5. The term 'periorbital' refers to the area around the eyes. Therefore, the physical finding directly points to a rash on the Eyelids.")
    
    # Evaluate Answer Choices
    print("\nEvaluating the choices:")
    print("A. Dorsum of the hands: Possible in Dermatomyositis (Gottron's papules), but not the finding described.")
    print("B. Nose: A rash here is more typical for lupus (malar rash).")
    print("C. Eyelids: This directly corresponds to the described 'periorbital erythema' and the classic Heliotrope rash.")
    print("D. Groin: Not a typical location for rashes associated with this systemic disease.")
    print("E. Shoulders: Possible in Dermatomyositis (Shawl sign), but not the finding described.")

    print("\nConclusion:")
    print("The most accurate answer is Eyelids, as this is the specific location of the described 'periorbital erythema' (Heliotrope rash).")

# Execute the analysis
solve_medical_case()