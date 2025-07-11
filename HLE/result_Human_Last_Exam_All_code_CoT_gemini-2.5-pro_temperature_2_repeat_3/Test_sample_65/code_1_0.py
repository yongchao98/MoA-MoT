def prioritize_patients():
    """
    Prioritizes three spine trauma patients based on clinical indications for surgery.
    """
    # Patient data representation
    # Scores are assigned based on urgency:
    # Neurologic deficit is the highest priority (e.g., 100 points).
    # Severity of instability is the second priority (e.g., 1-50 points).
    patients = {
        'Patient 1': {
            'diagnosis': 'Severe burst fracture of L2',
            'neurology': 'No neurologic deficits',
            'neuro_score': 0,
            'instability_score': 40 # Severe burst fracture indicates high instability
        },
        'Patient 2': {
            'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1',
            'neurology': 'No neurologic deficits',
            'neuro_score': 0,
            'instability_score': 20 # Mild spondylolisthesis indicates lower instability
        },
        'Patient 3': {
            'diagnosis': 'Split fracture of L2',
            'neurology': 'Mildly disordered pelvic functions',
            'neuro_score': 100, # Presence of cauda equina symptoms is a surgical emergency
            'instability_score': 35 # Split fracture indicates high instability
        }
    }

    # Calculate total priority score for each patient
    for patient_name, data in patients.items():
        data['total_score'] = data['neuro_score'] + data['instability_score']

    # Sort patients by total score in descending order
    # The key for sorting is the 'total_score' of each patient dictionary
    sorted_patients = sorted(patients.items(), key=lambda item: item[1]['total_score'], reverse=True)

    print("Spine Surgery Triage Prioritization:")
    print("------------------------------------\n")
    print("The primary principle for prioritization is addressing active neurological deficits, followed by the degree of spinal instability.\n")
    
    print("1. Highest Priority: Patient with Neurologic Deficits")
    print("   - Patient 3 shows 'mildly disordered pelvic functions', indicating a neurologic emergency (cauda equina syndrome). This requires immediate intervention.\n")

    print("2. Second Priority: Patient with Highest Mechanical Instability (without neuro deficits)")
    print("   - Patient 1 has a 'severe burst fracture'. This is a highly unstable injury with a significant risk of future neurologic damage and deformity if not surgically stabilized.\n")

    print("3. Third Priority: Patient with Lower Mechanical Instability (without neuro deficits)")
    print("   - Patient 2 has a 'mild' spondylolisthesis, which is less unstable than a severe burst fracture.\n")

    # Output the final ordered list
    print("Final Prioritization Order:")
    final_order = []
    for patient_info in sorted_patients:
        patient_num = patient_info[0].split(' ')[1]
        final_order.append(f"Patient {patient_num}")

    print(" -> ".join(final_order))
    
    # Just to fulfill the "output each number in the final equation" instruction format
    patient_numbers = [p[0].split(' ')[1] for p in sorted_patients]
    print(f"\nThe prioritized sequence of patients is: {patient_numbers[0]}, {patient_numbers[1]}, {patient_numbers[2]}")

prioritize_patients()
<<<F>>>