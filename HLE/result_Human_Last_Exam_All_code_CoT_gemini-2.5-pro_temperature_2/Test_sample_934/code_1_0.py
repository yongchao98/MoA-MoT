import collections

def analyze_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    
    # Patient Data from the vignette
    patient_data = {
        "age": 53,
        "symptoms": ["substernal chest pain (10/10)", "pain with swallowing"],
        "history": ["alcohol use disorder", "heavy smoker (2 packs/day for 20 years)"],
        "labs": ["elevated C-reactive protein", "leukocytosis"],
        "imaging": ["esophageal lumen narrowing", "wall thickening"],
        "endoscopy": "normal (no erythema, ulcers, plaques, or strictures)"
    }
    
    # Initialize scores for each diagnosis
    diagnosis_scores = collections.OrderedDict({
        "A. Streptococcal esophagitis": 0,
        "B. Esophageal adenocarcinoma": 0,
        "C. Esophageal squamous cell carcinoma": 0,
        "D. GERD": 0,
        "E. Herpes esophagitis": 0
    })
    
    print("Evaluating patient data against possible diagnoses:\n")
    
    # Scoring Logic
    
    # --- C. Esophageal squamous cell carcinoma (SCC) ---
    reasoning_scc = []
    if "heavy smoker (2 packs/day for 20 years)" in patient_data["history"]:
        diagnosis_scores["C. Esophageal squamous cell carcinoma"] += 3
        reasoning_scc.append("(+3) Major risk factor: Heavy smoking.")
    if "alcohol use disorder" in patient_data["history"]:
        diagnosis_scores["C. Esophageal squamous cell carcinoma"] += 3
        reasoning_scc.append("(+3) Major risk factor: Alcohol use disorder.")
    if "wall thickening" in patient_data["imaging"] and "esophageal lumen narrowing" in patient_data["imaging"]:
        diagnosis_scores["C. Esophageal squamous cell carcinoma"] += 2
        reasoning_scc.append("(+2) Imaging findings (wall thickening, narrowing) are classic for an infiltrative tumor.")
    if patient_data["endoscopy"] == "normal (no erythema, ulcers, plaques, or strictures)":
        diagnosis_scores["C. Esophageal squamous cell carcinoma"] += 2
        reasoning_scc.append("(+2) A normal endoscopy with positive imaging suggests a submucosal tumor, a known presentation of SCC.")
    print(f"Esophageal Squamous Cell Carcinoma Likelihood: {diagnosis_scores['C. Esophageal squamous cell carcinoma']}")
    for reason in reasoning_scc:
        print(f"  - {reason}")
    print("-" * 20)
    
    # --- A. Streptococcal esophagitis & E. Herpes esophagitis (Infectious) ---
    reasoning_infectious = []
    if "pain with swallowing" in patient_data["symptoms"]:
        diagnosis_scores["A. Streptococcal esophagitis"] += 1
        diagnosis_scores["E. Herpes esophagitis"] += 1
        reasoning_infectious.append("(+1) Symptom of odynophagia is present.")
    if patient_data["endoscopy"] == "normal (no erythema, ulcers, plaques, or strictures)":
        diagnosis_scores["A. Streptococcal esophagitis"] -= 3
        diagnosis_scores["E. Herpes esophagitis"] -= 3
        reasoning_infectious.append("(-3) Normal endoscopy makes infectious esophagitis (which typically shows ulcers/plaques) very unlikely.")
    print(f"Infectious Esophagitis (Strep/Herpes) Likelihood: {diagnosis_scores['A. Streptococcal esophagitis']}")
    for reason in reasoning_infectious:
        print(f"  - {reason}")
    print("-" * 20)
    
    # --- B. Esophageal adenocarcinoma ---
    reasoning_adeno = []
    if "GERD" not in patient_data["history"]:
        diagnosis_scores["B. Esophageal adenocarcinoma"] -= 1
        reasoning_adeno.append("(-1) No history of GERD, the main risk factor for adenocarcinoma.")
    if patient_data["endoscopy"] == "normal (no erythema, ulcers, plaques, or strictures)":
        diagnosis_scores["B. Esophageal adenocarcinoma"] -= 3
        reasoning_adeno.append("(-3) Normal endoscopy makes adenocarcinoma highly unlikely as it usually presents as a visible mass.")
    print(f"Esophageal Adenocarcinoma Likelihood: {diagnosis_scores['B. Esophageal adenocarcinoma']}")
    for reason in reasoning_adeno:
        print(f"  - {reason}")
    print("-" * 20)
    
    # --- D. GERD ---
    reasoning_gerd = []
    if "wall thickening" in patient_data["imaging"]:
        diagnosis_scores["D. GERD"] -= 2
        reasoning_gerd.append("(-2) Imaging findings are not typical for GERD.")
    if patient_data["endoscopy"] == "normal (no erythema, ulcers, plaques, or strictures)":
        diagnosis_scores["D. GERD"] -= 1
        reasoning_gerd.append("(-1) While non-erosive GERD exists, it would not explain the severe symptoms or imaging findings.")
    print(f"GERD Likelihood: {diagnosis_scores['D. GERD']}")
    for reason in reasoning_gerd:
        print(f"  - {reason}")
    print("-" * 20)
    
    # Determine the most likely diagnosis
    most_likely_diagnosis = max(diagnosis_scores, key=diagnosis_scores.get)
    
    print("\nConclusion:")
    print(f"The highest-scoring diagnosis is '{most_likely_diagnosis}' with a score of {diagnosis_scores[most_likely_diagnosis]}.")
    print("This is due to the strong correlation between the patient's major risk factors (smoking, alcohol) and the specific combination of imaging and endoscopic findings.")

analyze_clinical_case()