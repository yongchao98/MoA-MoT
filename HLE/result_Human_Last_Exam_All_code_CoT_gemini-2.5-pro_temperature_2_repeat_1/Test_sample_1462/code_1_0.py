def find_most_likely_defect():
    """
    This function models the diagnostic process for the given clinical case.
    It scores potential diagnoses based on their association with the patient's symptoms.
    """

    # 1. Define patient symptoms from the clinical vignette
    patient_symptoms = {
        "macrosomia": True,           # Weight of 12-lb 1-oz
        "respiratory_distress": True, # O2 sat 89% & fluid in lung
        "micrognathia": True,         # Small jaw
    }

    # 2. Create a knowledge base of disease-symptom associations.
    # Scores: 3 = Strong association, 2 = Moderate, 1 = Weak, 0 = No association
    disease_associations = {
        "A. Pulmonary hypoplasia": {
            "macrosomia": 0, "respiratory_distress": 3, "micrognathia": 0
        },
        "B. Pleuroperitoneal membrane defect": {
            "macrosomia": 0, "respiratory_distress": 3, "micrognathia": 1
        },
        "C. Ventral foregut budding defect": {
            "macrosomia": 0, "respiratory_distress": 2, "micrognathia": 1
        },
        "D. Maternal diabetes": {
            "macrosomia": 3, "respiratory_distress": 3, "micrognathia": 1
        },
        "E. Fistula": {
            "macrosomia": 0, "respiratory_distress": 2, "micrognathia": 1
        }
    }

    print("Calculating likelihood scores based on patient's symptoms...")
    print("Scoring Key: 3=Strong, 2=Moderate, 1=Weak, 0=None\n")

    best_match = ""
    max_score = -1

    # 3. Calculate the score for each potential diagnosis
    for disease, associations in disease_associations.items():
        current_score = 0
        equation_parts = []
        
        # Check only for the symptoms the patient has
        for symptom, is_present in patient_symptoms.items():
            if is_present:
                # Get the score for this symptom-disease link
                symptom_score = associations.get(symptom, 0)
                current_score += symptom_score
                # Add the number to our equation string
                equation_parts.append(str(symptom_score))
        
        # Build the final equation string for printing
        equation_str = " + ".join(equation_parts)
        
        print(f"Diagnosis: {disease}")
        # Print the simple "equation" used for scoring this diagnosis
        print(f"  Score Breakdown (Macrosomia + Resp. Distress + Micrognathia):")
        print(f"  Calculation: {equation_str} = {current_score}\n")

        # Keep track of the highest score
        if current_score > max_score:
            max_score = current_score
            best_match = disease

    print("---------------------------------------------------------")
    print(f"The most likely underlying cause is the one with the highest score.")
    print(f"Conclusion: {best_match}")
    print("---------------------------------------------------------")


# Run the diagnostic model
find_most_likely_defect()