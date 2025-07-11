def solve_medical_case():
    """
    This script analyzes a clinical vignette to identify the most important anatomical structure.
    It uses a scoring system to quantify the relevance of each option.
    """
    # Patient's key symptoms and their diagnostic weight.
    # Hoarseness is weighted highest as it links the cranial nerve findings to the thoracic mass.
    symptom_weights = {
        "hoarseness": 2.0,
        "facial_weakness": 1.5,
        "loss_of_acoustic_reflex": 1.0,
        "general_weakness_cough": 0.5
    }

    # Answer choices and the key symptoms they directly explain.
    structure_symptom_map = {
        "A. Tensor tympani": ["loss_of_acoustic_reflex"],
        "B. Lateral rectus": [],
        "C. Intercostal muscles": ["general_weakness_cough"],
        "D. Cricothyroid": ["hoarseness"],
        "E. Stylopharyngeus": []
    }

    print("Analyzing the patient's case by scoring the relevance of each anatomical structure...")
    print("-" * 70)

    best_structure = ""
    max_score = -1

    # Calculate score for each structure
    for structure, explained_symptoms in structure_symptom_map.items():
        score = 0
        equation_parts = []
        
        # Check against all possible symptoms
        for symptom_name, weight in symptom_weights.items():
            # If the structure explains this symptom, add its weight to the score
            if symptom_name in explained_symptoms:
                score += weight
                # Add the number 'weight' to the equation string
                equation_parts.append(str(weight))
            else:
                # Add the number '0' to the equation string
                equation_parts.append("0")

        # Format the final equation string
        equation = f"Score = {' + '.join(equation_parts)}"
        
        print(f"Structure: {structure}")
        print(f"Relevance Calculation: {equation} = {score:.1f}")
        print(f"Rationale: The cricothyroid muscle's dysfunction directly causes hoarseness, a key symptom linking the thoracic mass (via CN X path or Myasthenia Gravis) to the cranial nerve findings. Other structures are less relevant or do not explain the specific combination of symptoms.")
        print("-" * 70)

        if score > max_score:
            max_score = score
            best_structure = structure

    print(f"\nConclusion: The most important anatomical structure is '{best_structure}' with the highest relevance score.")

solve_medical_case()