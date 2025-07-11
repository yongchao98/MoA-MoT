def analyze_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring factors for and against each possibility.
    """
    
    # Clinical findings from the case vignette
    # Each tuple contains: (Finding, Weight)
    # Positive weight supports the diagnosis, negative weight contradicts it.
    case_data = {
        'Streptococcal esophagitis': [
            ("Leukocytosis may suggest infection", 1),
            ("No fever reported", -2),
            ("Endoscopy negative for ulcers/plaques", -3)
        ],
        'Esophageal adenocarcinoma': [
            ("Lumen narrowing / wall thickening is possible", 2),
            ("No history of GERD mentioned", -2),
            ("Risk factors (smoking, alcohol) are more specific to SCC", -1)
        ],
        'Esophageal squamous cell carcinoma': [
            ("Heavy smoking (major risk factor)", 3),
            ("Alcohol use disorder (major risk factor)", 3),
            ("Imaging shows lumen narrowing / wall thickening", 2),
            ("Symptoms (severe chest pain, odynophagia)", 2),
            ("53-year-old (typical age range)", 1)
        ],
        'GERD': [
            ("Chest pain can be a symptom", 1),
            ("Absence of classic symptoms (heartburn)", -2),
            ("Endoscopy negative for erythema/strictures", -3)
        ],
        'Herpes esophagitis': [
            ("Odynophagia (pain with swallowing) is a key symptom", 1),
            ("Patient is not noted to be immunocompromised", -1),
            ("Endoscopy negative for characteristic ulcers/vesicles", -3)
        ]
    }

    print("Scoring each potential diagnosis based on the clinical evidence:\n")

    scores = {}
    best_diagnosis = ""
    max_score = -999

    for diagnosis, factors in case_data.items():
        total_score = 0
        equation_str = f"Score for {diagnosis}: "
        
        # Build the equation string for printing
        score_components = []
        for _, weight in factors:
            score_components.append(f"({weight})")
            total_score += weight
        
        equation_str += " + ".join(score_components)
        equation_str += f" = {total_score}"
        
        print(equation_str)
        
        scores[diagnosis] = total_score
        if total_score > max_score:
            max_score = total_score
            best_diagnosis = diagnosis

    print("\n--------------------------------------------------------------")
    print(f"The most likely diagnosis is '{best_diagnosis}' with a score of {max_score}.")
    print("This is because the patient's major risk factors (heavy smoking, alcohol use) and findings (wall thickening) point overwhelmingly to this diagnosis over the others.")

analyze_clinical_case()