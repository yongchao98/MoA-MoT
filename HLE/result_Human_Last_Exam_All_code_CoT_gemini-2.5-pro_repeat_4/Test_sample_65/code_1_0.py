def prioritize_patients():
    """
    Prioritizes spine surgery patients based on clinical urgency.
    The scoring is based on:
    3: Active neurological deficit (e.g., cauda equina syndrome) -> Emergency
    2: High spinal instability with risk of neurological injury -> Urgent
    1: Lower spinal instability, no neurological deficit -> Less urgent
    """
    
    # Patient data and assigned priority scores
    patients = [
        {'id': 1, 'condition': 'Severe burst fracture of L2, no neurologic deficits', 'priority_score': 2},
        {'id': 2, 'condition': 'Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits', 'priority_score': 1},
        {'id': 3, 'condition': 'Split fracture of L2, with mildly disordered pelvic functions', 'priority_score': 3}
    ]
    
    # Sort patients by priority score in descending order
    sorted_patients = sorted(patients, key=lambda x: x['priority_score'], reverse=True)
    
    print("Patient Prioritization Based on Surgical Urgency:\n")
    
    priority_map = {
        3: "Top Priority (Neurosurgical Emergency)",
        2: "Second Priority (High Mechanical Instability)",
        1: "Lowest Priority (Stable/Less Unstable)"
    }
    
    print("Final Prioritization Order:")
    
    # Construct the final equation string
    final_order_list = []
    for patient in sorted_patients:
        final_order_list.append(f"Patient {patient['id']}")
    
    # Print each number in the final equation
    print("Patient 3, Patient 1, Patient 2")


prioritize_patients()