def diagnose_and_find_test():
    """
    This script simulates the clinical reasoning process to determine the most likely diagnosis
    and the corresponding confirmatory physical exam test.
    """

    # 1. Define patient's key clinical features
    patient_symptoms = {
        "pain_location": "L4-S1 distribution",
        "aggravating_factor": "pain worse when supine",
        "exam_position": "left decubitus",
        "test_type": "resisted movement"
    }

    print("Step 1: Analyzing patient's clinical presentation.")
    print(f"- Symptoms: {patient_symptoms['pain_location']}, {patient_symptoms['aggravating_factor']}")
    print(f"- Exam Setup: Patient in {patient_symptoms['exam_position']}, test involves {patient_symptoms['test_type']}\n")

    # 2. Define a knowledge base of possible diagnoses and their associated tests
    # Note: Actions are described for the right leg, with the patient in the left decubitus position.
    knowledge_base = {
        "Piriformis Syndrome": {
            "symptoms": ["L4-S1 distribution", "pain worse when supine"],
            "confirmatory_test": "External Rotation",
            "test_type": "resisted movement",
            "rationale": "The piriformis is a primary external rotator. Resisted external rotation causes it to contract, compressing the sciatic nerve."
        },
        "SI Joint Dysfunction": {
            "symptoms": ["L4-S1 distribution", "pain worse when supine"],
            "confirmatory_test": "Extension",
            "test_type": "passive movement", # Primarily passive (Gaenslen's test), but resisted extension can stress it.
            "rationale": "Tests like Gaenslen's involve passive hyperextension to stress the SI joint."
        },
        "Gluteal Tendinopathy": {
            "symptoms": ["L4-S1 distribution"], # Can mimic an L5 pattern
            "confirmatory_test": "Abduction",
            "test_type": "resisted movement",
            "rationale": "Resisted abduction stresses the gluteus medius/minimus tendons."
        },
        "Lumbar Disc Herniation": {
            "symptoms": ["L4-S1 distribution"],
            "confirmatory_test": "Flexion", # As in a Straight Leg Raise, but less specific here.
            "test_type": "passive movement",
            "rationale": "Typically diagnosed with tests like the Straight Leg Raise (hip flexion) in a supine position."
        }
    }
    
    print("Step 2: Evaluating potential diagnoses based on symptoms.")
    scores = {}
    equations = {}
    for diagnosis, data in knowledge_base.items():
        score = 0
        equation_str = []
        # Score based on matching symptoms
        if patient_symptoms["pain_location"] in data["symptoms"]:
            score += 1
            equation_str.append("1 (for L4-S1 pain)")
        if patient_symptoms["aggravating_factor"] in data["symptoms"]:
            score += 1
            equation_str.append("1 (for pain worse when supine)")
        
        # Bonus point if the test type matches the scenario
        if patient_symptoms["test_type"] == data["test_type"]:
            score += 1
            equation_str.append(f"1 (for matching test type: {data['test_type']})")

        scores[diagnosis] = score
        equations[diagnosis] = " + ".join(equation_str)
    
    # Print scoring
    for diagnosis, score in scores.items():
        print(f"- {diagnosis} Score Equation: {equations[diagnosis]} = {score}")

    # 3. Determine the most likely diagnosis
    most_likely_diagnosis = max(scores, key=scores.get)
    print(f"\nStep 3: Identifying the most likely diagnosis.")
    print(f"The diagnosis with the highest score is: {most_likely_diagnosis} (Score: {scores[most_likely_diagnosis]})\n")

    # 4. Identify the confirmatory test for the top diagnosis
    final_test = knowledge_base[most_likely_diagnosis]['confirmatory_test']
    rationale = knowledge_base[most_likely_diagnosis]['rationale']
    answer_choices = {
        'A': 'Abduction',
        'B': 'Adduction',
        'C': 'Internal Rotation',
        'D': 'External Rotation',
        'E': 'Flexion',
        'F': 'Extension'
    }

    final_answer_letter = [letter for letter, action in answer_choices.items() if action == final_test][0]

    print("Step 4: Determining the confirmatory action.")
    print(f"The confirmatory test for {most_likely_diagnosis} is resisted '{final_test}'.")
    print(f"Rationale: {rationale}")
    print(f"\nThis corresponds to answer choice {final_answer_letter}.")


if __name__ == '__main__':
    diagnose_and_find_test()