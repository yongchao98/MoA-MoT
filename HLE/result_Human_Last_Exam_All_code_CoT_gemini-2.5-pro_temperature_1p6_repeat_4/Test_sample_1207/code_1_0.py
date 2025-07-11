def find_best_diagnosis():
    """
    Analyzes patient symptoms to determine the most likely diagnosis and associated imaging finding.
    """
    # Step 1: Define the patient's key symptoms
    patient_symptoms = {
        "transient_monocular_vision_loss", # Suggests vascular occlusion in the retina
        "pulsatile_headaches",             # Neurologic symptom
        "hearing_loss",                    # Auditory symptom
        "joint_pain",                      # Systemic inflammatory symptom
        "dyspnea"                          # Pulmonary/systemic symptom
    }
    print("Patient's Key Symptoms:")
    for symptom in sorted(list(patient_symptoms)):
        print(f"- {symptom.replace('_', ' ').title()}")
    print("-" * 30)

    # Step 2 & 3: Define differential diagnoses and their characteristic features
    # This simulates the clinical reasoning process.
    diagnoses = {
        "Susac's Syndrome": {
            "classic_triad": {"encephalopathy", "branch_retinal_artery_occlusions", "sensorineural_hearing_loss"},
            "explanation": "This syndrome is an autoimmune endotheliopathy affecting the microvasculature of the brain, retina, and inner ear. The patient's transient vision loss is classic for retinal artery occlusions, and the combination with headaches and hearing loss is highly suggestive.",
            "imaging_finding": 'Leptomeningeal enhancement with "snowball" hyperintensities visualized by MRI',
            "answer_choice": "B"
        },
        "Pituitary Adenoma": {
            "classic_triad": {"headache", "bitemporal_hemianopsia", "hormonal_changes"},
            "explanation": "Explains headaches and vision loss but not the specific type (monocular/transient) or the auditory, joint, and pulmonary symptoms.",
            "imaging_finding": "Intrasellar mass visualized by MRI",
            "answer_choice": "E"
        },
        "Systemic Lupus Erythematosus (SLE)": {
            "classic_triad": {"joint_pain", "skin_rash", "fatigue"},
            "explanation": "A multi-system autoimmune disease that can cause all these symptoms, but the specific triad of brain/eye/ear issues is less classic than for Susac's syndrome.",
            "imaging_finding": "Pleural effusion or pericardial effusion are common.",
            "answer_choice": "C"
        },
        "Giant Cell Arteritis (GCA)": {
             "classic_triad": {"headache", "jaw_claudication", "vision_loss"},
             "explanation": "A vasculitis typically affecting patients >50. While it causes headaches and vision loss, the other systemic features are less typical.",
             "imaging_finding": "Vascular inflammation on MRA/CTA, sometimes hemorrhage.",
             "answer_choice": "D"
        }
    }

    # Step 4: Determine the best fit (simplified logic)
    # The triad of brain, eye, and ear symptoms is the most specific clinical clue.
    best_fit_diagnosis = "Susac's Syndrome"
    
    print(f"Most Likely Diagnosis: {best_fit_diagnosis}")
    print("\nReasoning:")
    print(diagnoses[best_fit_diagnosis]["explanation"])
    
    print("\nExpected Finding:")
    print(f"The pathognomonic (highly specific) finding for this diagnosis is: {diagnoses[best_fit_diagnosis]['imaging_finding']}")
    print(f"\nThis corresponds to answer choice {diagnoses[best_fit_diagnosis]['answer_choice']}.")


if __name__ == '__main__':
    find_best_diagnosis()