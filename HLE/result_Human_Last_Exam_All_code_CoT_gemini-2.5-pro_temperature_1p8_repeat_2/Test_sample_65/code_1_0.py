import sys

def triage_spine_patients():
    """
    This script prioritizes three spine surgery patients based on their diagnosis
    and the urgency of their surgical indications.
    """
    # Step 1: Define patient data
    # Each patient is represented as a dictionary containing their ID, diagnosis,
    # and key clinical factors for scoring.
    patients = [
        {'id': 1, 'description': 'Severe burst fracture of L2, no neurologic deficits.', 'has_neuro_deficit': False, 'instability_level': 8},
        {'id': 2, 'description': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.', 'has_neuro_deficit': False, 'instability_level': 4},
        {'id': 3, 'description': 'Split fracture of L2, with mildly disordered pelvic functions.', 'has_neuro_deficit': True, 'instability_level': 6}
    ]

    # Step 2: Explain the prioritization logic
    print("Spine Surgery Triage Plan:")
    print("="*30)
    print("Patients are prioritized based on two main factors:")
    print("1. Presence of Neurologic Deficits: This is the most critical factor and receives the highest score.")
    print("2. Degree of Spinal Instability: Severe instability poses a high risk of future injury.")
    print("="*30 + "\n")

    # Step 3: Calculate a priority score for each patient
    # The scoring reflects the clinical urgency.
    for patient in patients:
        neuro_score = 100 if patient['has_neuro_deficit'] else 0
        instability_score = patient['instability_level']
        patient['priority_score'] = neuro_score + instability_score

    # Step 4: Sort patients by their priority score in descending order
    sorted_patients = sorted(patients, key=lambda p: p['priority_score'], reverse=True)

    # Step 5: Print the detailed results and the final priority equation
    print("Triage Results (from highest to lowest priority):\n")
    
    patient_order = []
    for i, patient in enumerate(sorted_patients):
        print(f"Priority #{i+1}: Patient {patient['id']}")
        print(f"  - Diagnosis: {patient['description']}")
        print(f"  - Rationale: Score={patient['priority_score']} (Neuro Deficit: {patient['has_neuro_deficit']}, Instability: {patient['instability_level']})")
        print("-" * 25)
        patient_order.append(str(patient['id']))

    # Fulfills the requirement to output each number in a final equation format
    print("\nFinal Prioritization Order Equation:")
    print(f"Patient {patient_order[0]} > Patient {patient_order[1]} > Patient {patient_order[2]}")


# Run the triage function
triage_spine_patients()

# Based on the logic, the correct order is Patient 3, Patient 1, Patient 2.
# This corresponds to answer choice F.
sys.stdout.write("<<<F>>>")