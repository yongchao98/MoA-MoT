def prioritize_spine_patients():
    """
    Prioritizes three spine trauma patients based on clinical urgency.
    """
    # Patient data with diagnosis and neurologic status
    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture', 'neurology': 'no deficits'},
        {'id': 2, 'diagnosis': 'Compression fracture with mild spondylolisthesis', 'neurology': 'no deficits'},
        {'id': 3, 'diagnosis': 'Split fracture', 'neurology': 'mildly disordered pelvic functions'}
    ]

    # Assign priority scores based on clinical hierarchy:
    # 3 = Active Neurologic Deficit (Surgical Emergency)
    # 2 = Severe Mechanical Instability (Urgent)
    # 1 = Mechanical Instability (Less Urgent)
    for patient in patients:
        if patient['neurology'] != 'no deficits':
            patient['priority_score'] = 3
        elif 'Severe burst fracture' in patient['diagnosis']:
            patient['priority_score'] = 2
        else:
            patient['priority_score'] = 1

    # Sort patients in descending order of priority
    sorted_patients = sorted(patients, key=lambda p: p['priority_score'], reverse=True)

    print("Prioritization of Patients for Surgery (Top to Lowest Priority):")
    
    # Extracting the final order of patient numbers for the "equation"
    final_order_numbers = [p['id'] for p in sorted_patients]
    
    # Print the detailed breakdown
    for i, patient in enumerate(sorted_patients):
        print(f"Priority {i+1}: Patient {patient['id']} (Condition: {patient['diagnosis']}, Status: {patient['neurology']})")
    
    # Print the final ordered equation as requested
    print("\nFinal Priority Order:")
    # Using ' -> ' to represent the flow of priority
    equation_output = f"Patient {final_order_numbers[0]} -> Patient {final_order_numbers[1]} -> Patient {final_order_numbers[2]}"
    print(equation_output)

prioritize_spine_patients()
<<<F>>>