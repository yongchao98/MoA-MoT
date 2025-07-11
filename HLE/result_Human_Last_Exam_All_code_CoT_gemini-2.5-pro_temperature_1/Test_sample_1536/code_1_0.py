def solve_clinical_case():
    """
    This script analyzes a clinical case to identify the most relevant anatomical structure.
    It scores each potential answer based on how well it explains the patient's symptoms.
    """
    symptoms = {
        "facial_weakness_CN_VII": 2,
        "loss_of_acoustic_reflex_CN_VII": 2,
        "hoarseness_cough_CN_X_larynx": 3,
        "thoracic_mass_impact_on_recurrent_laryngeal_nerve": 3
    }

    # Anatomical structures and their primary nerve/function
    # We will score each based on the patient's symptoms
    choices = {
        "A": {"name": "Tensor tympani", "nerve": "CN V", "function": "Tenses eardrum"},
        "B": {"name": "Lateral rectus", "nerve": "CN VI", "function": "Eye abduction"},
        "C": {"name": "Intercostal muscles", "nerve": "Intercostal nerves", "function": "Respiration"},
        "D": {"name": "Cricothyroid", "nerve": "CN X (Superior Laryngeal)", "function": "Laryngeal muscle - voice/pitch"},
        "E": {"name": "Stylopharyngeus", "nerve": "CN IX", "function": "Pharynx elevation"}
    }

    scores = {key: 0 for key in choices}
    reasoning = {key: [] for key in choices}

    # Scoring Logic
    # A. Tensor tympani
    # No relation to CN VII, CN X, or laryngeal symptoms.
    scores["A"] = 0
    reasoning["A"].append("Not related to Facial (CN VII) or Vagus (CN X) nerve symptoms.")

    # B. Lateral rectus
    # No relation to symptoms.
    scores["B"] = 0
    reasoning["B"].append("Patient has no reported eye movement issues (CN VI).")

    # C. Intercostal muscles
    # Located in thorax, but doesn't explain cranial nerve signs.
    scores["C"] = 1
    reasoning["C"].append("Located in the thorax but does not explain hoarseness or facial weakness.")

    # D. Cricothyroid
    # Directly related to hoarseness. Part of the laryngeal system innervated by the Vagus nerve (CN X),
    # whose recurrent laryngeal branch is implicated by the thoracic mass.
    scores["D"] += symptoms["hoarseness_cough_CN_X_larynx"]
    reasoning["D"].append(f"Directly related to hoarseness (laryngeal muscle). Score +{symptoms['hoarseness_cough_CN_X_larynx']}.")
    scores["D"] += symptoms["thoracic_mass_impact_on_recurrent_laryngeal_nerve"]
    reasoning["D"].append(f"The Vagus nerve system, which innervates this muscle, is affected by the thoracic mass. Score +{symptoms['thoracic_mass_impact_on_recurrent_laryngeal_nerve']}.")

    # E. Stylopharyngeus
    # Minor relation to throat function, but not the primary muscle for voice.
    scores["E"] = 0
    reasoning["E"].append("Innervated by CN IX; less relevant to the primary symptoms than CN X.")

    print("Analyzing the patient's symptoms to identify the most important anatomical structure...\n")
    for key, choice in choices.items():
        print(f"Choice {key}: {choice['name']}")
        print(f"  - Function: {choice['function']}")
        print(f"  - Reasoning: {' '.join(reasoning[key])}")
        print(f"  - Final Score: {scores[key]}\n")

    # Find the best choice
    best_choice_key = max(scores, key=scores.get)
    best_choice_name = choices[best_choice_key]['name']

    print("--- Conclusion ---")
    print(f"The structure that best explains the key symptom of hoarseness, which is linked to the thoracic mass via the Vagus nerve system, is the {best_choice_name}.")
    print(f"Therefore, the most important anatomical structure among the choices is '{best_choice_name}'.")

solve_clinical_case()
<<<D>>>