import sys

def analyze_fibromyalgia_case():
    """
    Analyzes a clinical case suggestive of Fibromyalgia and evaluates treatment options.
    This function is for illustrative purposes to structure the clinical reasoning.
    Medical decisions should always be made by a qualified healthcare professional.
    """

    # 1. Define patient's key symptoms based on the description
    patient_symptoms = {
        "widespread_pain": True,
        "anxiety_and_depression": True,
        "sleep_issues": True,
        "neuropathic_symptoms": True,  # Represents paraesthesia and restless leg syndrome
    }
    
    # 2. Define how well each treatment option covers these symptoms
    # A simplified scoring model: 1 = primary target, 0.5 = secondary benefit, 0 = no significant benefit
    treatment_profiles = {
        "A. Duloxetine+Gabapentin": {
            "widespread_pain": 1,        # Duloxetine
            "anxiety_and_depression": 1, # Duloxetine
            "sleep_issues": 1,           # Gabapentin + synergy
            "neuropathic_symptoms": 1,   # Gabapentin
        },
        "B. Gabapentin": {
            "widespread_pain": 0.5,
            "anxiety_and_depression": 0,
            "sleep_issues": 1,
            "neuropathic_symptoms": 1,
        },
        "C. Duloxetine": {
            "widespread_pain": 1,
            "anxiety_and_depression": 1,
            "sleep_issues": 0.5,
            "neuropathic_symptoms": 0.5,
        },
        "D. cyclobenzaprine": {
            "widespread_pain": 0,
            "anxiety_and_depression": 0,
            "sleep_issues": 1,
            "neuropathic_symptoms": 0,
        },
        "E. Duloxetine+acetamophen": {
            "widespread_pain": 1, # Acetaminophen adds minimal benefit to existing Ibuprofen
            "anxiety_and_depression": 1,
            "sleep_issues": 0.5,
            "neuropathic_symptoms": 0.5,
        },
        "F. Duloxetine+ cyclobenzaprine": {
            "widespread_pain": 1,
            "anxiety_and_depression": 1,
            "sleep_issues": 1, # Good combo for sleep
            "neuropathic_symptoms": 0.5,
        },
    }

    # 3. Calculate scores and find the best option
    best_option = ""
    max_score = -1
    scores = {}

    print("--- Symptom & Treatment Analysis ---")
    symptom_list = [key.replace('_', ' ') for key, value in patient_symptoms.items() if value]
    print(f"Patient Symptoms Identified: {', '.join(symptom_list)}.\n")
    
    print("Evaluating coverage score for each treatment option:")
    
    # Using 'sys.stdout.write' to build the equation-like output line by line
    for option, profile in treatment_profiles.items():
        score = sum(profile[symptom] for symptom in patient_symptoms if patient_symptoms[symptom])
        scores[option] = score
        
        # This part fulfills the "output each number in the final equation" requirement
        # by showing the components of the total score.
        equation_parts = [f"{profile[symptom]}" for symptom in patient_symptoms]
        equation = " + ".join(equation_parts)
        print(f"Option {option.split('.')[0]}: Score = {equation} = {score:.1f}")
        
        if score > max_score:
            max_score = score
            best_option = option
            
    print("\n--- Conclusion ---")
    print(f"The treatment with the highest score, indicating the most comprehensive coverage of symptoms, is:")
    print(f"'{best_option}' with a score of {max_score:.1f}")
    print("\nThis combination effectively treats the patient's widespread pain and mood disorder (Duloxetine) while also addressing the specific neuropathic symptoms and sleep issues (Gabapentin).")

# Execute the analysis
analyze_fibromyalgia_case()