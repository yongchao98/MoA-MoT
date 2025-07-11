def diagnose_esophageal_condition():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    This function uses a scoring system based on how well the patient's data
    matches the characteristics of each possible diagnosis.
    """

    # Patient Data Points:
    # (+) indicates a supporting feature, (-) indicates a contradictory feature
    # Patient has: Heavy smoking, alcohol use, chest pain, odynophagia,
    # imaging shows narrowing/thickening, BUT endoscopy shows no ulcers/plaques/erythema.

    diagnoses = {
        "A": {"name": "Streptococcal esophagitis", "score": 0, "reasoning": []},
        "B": {"name": "Esophageal adenocarcinoma", "score": 0, "reasoning": []},
        "C": {"name": "Esophageal squamous cell carcinoma", "score": 0, "reasoning": []},
        "D": {"name": "GERD", "score": 0, "reasoning": []},
        "E": {"name": "Herpes esophagitis", "score": 0, "reasoning": []}
    }

    # --- Scoring Logic ---

    # A. Streptococcal esophagitis: Rare, typically causes visible exudates.
    score_a = -2  # (-2 for normal endoscopy)
    diagnoses["A"]["score"] = score_a
    diagnoses["A"]["reasoning"].append(f"Normal endoscopy contradicts visible exudates/ulcers: {score_a}")

    # B. Esophageal adenocarcinoma: Linked to GERD (not reported), usually presents as a visible mass.
    score_b = 1 - 1 # (+1 for smoking as a weak risk factor, -1 for normal endoscopy being atypical)
    diagnoses["B"]["score"] = score_b
    diagnoses["B"]["reasoning"].append(f"Smoking is a weak risk factor (+1), but normal endoscopy is atypical (-1). Total: {score_b}")


    # C. Esophageal squamous cell carcinoma (SCC): Strong link to smoking and alcohol. Can grow submucosally.
    score_c = 3 + 3 + 3  # (+3 for major risk factor: smoking, +3 for major risk factor: alcohol, +3 for matching the unique presentation)
    diagnoses["C"]["score"] = score_c
    diagnoses["C"]["reasoning"].append(f"Heavy smoking risk factor (+3) + Alcohol use risk factor (+3) + Fits pattern of submucosal growth (imaging changes with normal endoscopy) (+3). Total: {score_c}")

    # D. GERD: Unlikely with severe symptoms and a normal endoscopy.
    score_d = 1 - 3 # (+1 for chest pain, -3 for normal endoscopy with severe symptoms)
    diagnoses["D"]["score"] = score_d
    diagnoses["D"]["reasoning"].append(f"Can cause chest pain (+1), but normal endoscopy with 10/10 pain is highly unlikely (-3). Total: {score_d}")

    # E. Herpes esophagitis: Typically presents with "punched-out" ulcers on endoscopy.
    score_e = -3 # (-3 for normal endoscopy)
    diagnoses["E"]["score"] = score_e
    diagnoses["E"]["reasoning"].append(f"Normal endoscopy contradicts typical 'punched-out' ulcers: {score_e}")


    # --- Print Results ---
    print("Evaluating diagnoses based on a scoring system:")
    highest_score = -float('inf')
    most_likely_diagnosis_key = None

    for key, data in diagnoses.items():
        print(f"\nDiagnosis ({key}) {data['name']}:")
        print(" | ".join(data['reasoning']))
        print(f"Final Score: {data['score']}")
        if data['score'] > highest_score:
            highest_score = data['score']
            most_likely_diagnosis_key = key

    print("\n-------------------------------------------------")
    print(f"The most likely diagnosis is ({most_likely_diagnosis_key}) {diagnoses[most_likely_diagnosis_key]['name']} with a score of {highest_score}.")
    print("This is because the patient's major risk factors (smoking, alcohol) and the specific clinical finding of a submucosal mass (wall thickening on imaging with a normal endoscopic appearance) are classic for this condition.")

    # This is the final answer for the user's prompt.
    print("<<<C>>>")

if __name__ == '__main__':
    diagnose_esophageal_condition()