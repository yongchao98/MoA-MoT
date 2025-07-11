import collections

def analyze_clinical_case():
    """
    Analyzes a clinical case by scoring potential diagnoses based on textual evidence.
    This is a logical analysis tool, not a substitute for professional medical advice.
    """
    # Step 1: Define the clinical information from the case study
    patient_data = {
        "timeline_days": 29,
        "procedure": "Whipple procedure", # A major abdominal surgery
        "symptoms": ["oxygen level 82% on 3L", "crackles on both sides", "gasping for air"],
        "history": ["major blood loss", "blood transfusions"]
    }

    # Step 2: Define the possible causes (answer choices)
    options = {
        'A': 'Acute blood transfusion reaction',
        'B': 'Iodine-related reaction',
        'C': 'Sensitivity reaction',
        'D': 'Sepsis',
        'E': 'Myocyte necrosis',
        'F': 'Respiratory deconditioning',
        'G': 'Lung exhaustion',
        'H': 'Air pollution sensitivity'
    }

    # Step 3: Score each option based on a simple ruleset
    scores = collections.defaultdict(int)
    reasoning = collections.defaultdict(list)

    # Key clinical features derived from the patient data
    is_major_surgery = "Whipple procedure" in patient_data["procedure"]
    is_long_timeline = patient_data["timeline_days"] > 7
    is_acute_timeline = patient_data["timeline_days"] <= 1
    has_transfusion_history = "blood transfusions" in patient_data["history"]
    has_ards_symptoms = "crackles on both sides" in patient_data["symptoms"] and "oxygen level 82%" in patient_data["symptoms"][0]

    # --- Scoring Logic ---

    # Option A: Acute blood transfusion reaction
    if has_transfusion_history and has_ards_symptoms:
        scores['A'] += 3
        reasoning['A'].append("+3 for transfusion history and ARDS-like symptoms (could be TRALI/TACO).")
    if not is_acute_timeline:
        scores['A'] -= 4
        reasoning['A'].append(f"-4 because the term 'Acute' is a poor fit for a {patient_data['timeline_days']}-day timeline.")

    # Option B: Iodine-related reaction
    if not is_acute_timeline:
        scores['B'] -= 3
        reasoning['B'].append(f"-3 because a reaction would be immediate, not {patient_data['timeline_days']} days later.")
    
    # Option C: Sensitivity reaction
    scores['C'] -= 2
    reasoning['C'].append("-2 for being a non-specific diagnosis.")

    # Option D: Sepsis
    if is_major_surgery:
        scores['D'] += 3
        reasoning['D'].append("+3 as major surgery is a high risk factor for post-op infection.")
    if is_long_timeline:
        scores['D'] += 2
        reasoning['D'].append(f"+2 as the {patient_data['timeline_days']}-day timeline fits the development of a serious infection.")
    if has_ards_symptoms:
        scores['D'] += 3
        reasoning['D'].append("+3 as sepsis is a primary cause of ARDS, which matches the lung symptoms.")

    # Option E: Myocyte necrosis (Heart Attack)
    if has_ards_symptoms:
        scores['E'] += 1
        reasoning['E'].append("+1 as a heart attack can cause pulmonary edema, but the primary signs are respiratory.")

    # Options F, G, H: Deconditioning, Exhaustion, Pollution
    scores['F'] -= 3
    reasoning['F'].append("-3 as deconditioning does not cause such severe acute hypoxemia and bilateral crackles.")
    scores['G'] -= 4
    reasoning['G'].append("-4 as 'Lung exhaustion' is not a standard medical term.")
    scores['H'] -= 4
    reasoning['H'].append("-4 as this is an unlikely cause of this acute presentation in a post-op patient.")

    # Step 4: Present the final analysis
    print("Clinical Case Analysis (Logical Scoring based on Text):\n")
    
    # Sort options by score for a ranked list
    sorted_options = sorted(options.keys(), key=lambda k: scores[k], reverse=True)

    for key in sorted_options:
        print(f"Option {key}: {options[key]}")
        print(f"  Score: {scores[key]}")
        print(f"  Reasoning: {' '.join(reasoning[key])}")
        print("-" * 30)

if __name__ == '__main__':
    analyze_clinical_case()