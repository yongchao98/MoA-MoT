import pandas as pd

def prioritize_patients():
    """
    This function prioritizes three spine surgery patients based on their clinical presentation.
    It uses a scoring system to rank them according to surgical urgency.
    """
    
    # Patient data
    patients_data = {
        'Patient 1': {
            'diagnosis': 'Severe burst fracture of L2',
            'neuro_status': 'No neurologic deficits',
            'id': 1
        },
        'Patient 2': {
            'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1',
            'neuro_status': 'No neurologic deficits',
            'id': 2
        },
        'Patient 3': {
            'diagnosis': 'Split fracture of L2',
            'neuro_status': 'Mildly disordered pelvic functions',
            'id': 3
        }
    }

    # Prioritization logic and scoring
    # A higher score indicates higher priority.
    # Neurologic deficit (especially cauda equina signs) is the most critical factor.
    # Fracture instability is the second most critical factor.
    scores = {}
    reasoning = {}
    
    for name, data in patients_data.items():
        score = 0
        reasons = []
        
        # Rule 1: Neurologic Status
        if 'disordered pelvic functions' in data['neuro_status']:
            score += 100
            reasons.append("Urgent priority due to signs of cauda equina syndrome (disordered pelvic functions), a surgical emergency.")
        
        # Rule 2: Fracture Instability
        if 'Severe burst fracture' in data['diagnosis']:
            score += 50
            reasons.append("High priority due to a highly unstable severe burst fracture, which risks future neurologic injury and deformity.")
        elif 'Split fracture' in data['diagnosis']:
            # Instability is present, but neuro status is the primary driver for this patient
            score += 40
            reasons.append("Moderate priority based on fracture pattern (instability), but neuro status is the main concern.")
        elif 'mild traumatic spondylolisthesis' in data['diagnosis']:
            score += 20
            reasons.append("Lower priority due to a more stable fracture pattern (compression Fx with only mild listhesis) and no neuro deficits.")
            
        scores[name] = score
        reasoning[name] = " ".join(reasons)

    # Sort patients by score in descending order
    sorted_patients = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    
    # Print the prioritized list and the reasoning
    print("Surgical Prioritization Order:\n")
    priority = 1
    final_order = []
    for patient_name, score in sorted_patients:
        patient_id = patients_data[patient_name]['id']
        final_order.append(str(patient_id))
        print(f"Priority #{priority}: Patient {patient_id}")
        print(f"Diagnosis: {patients_data[patient_name]['diagnosis']}, {patients_data[patient_name]['neuro_status']}")
        print(f"Priority Score: {score}")
        print(f"Reasoning: {reasoning[patient_name]}\n")
        priority += 1
        
    final_equation = " -> ".join(final_order)
    print("--------------------------------------------------")
    print(f"Final Priority Order (by Patient ID): {final_equation}")
    print("This corresponds to Answer Choice F.")
    print("--------------------------------------------------")


prioritize_patients()
<<<F>>>