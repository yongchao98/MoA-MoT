def prioritize_spine_patients():
    """
    Prioritizes spine surgery patients based on diagnosis and neurologic status.
    
    The prioritization logic is as follows:
    1.  Primary key: Presence of a neurologic deficit (True is higher priority).
    2.  Secondary key: Degree of mechanical instability (higher score is higher priority).
    """
    
    # Define patients with their conditions and priority scores
    # neuro_deficit: True if present, False if absent
    # instability_score: A numerical representation of fracture instability (higher is more unstable)
    # 3: Severe Burst Fracture (highly unstable)
    # 2: Split Fracture (a type of unstable burst fracture)
    # 1: Compression Fx with Mild Spondylolisthesis (moderately unstable)
    
    patients = [
        {'id': 1, 'diagnosis': 'Severe burst fracture of L2', 'neuro_deficit': False, 'instability_score': 3},
        {'id': 2, 'diagnosis': 'Compression fracture of L2 with mild traumatic spondylolisthesis', 'neuro_deficit': False, 'instability_score': 1},
        {'id': 3, 'diagnosis': 'Split fracture of L2, with mildly disordered pelvic functions', 'neuro_deficit': True, 'instability_score': 2}
    ]

    # Sort patients: highest priority first.
    # Python's sort is stable. We sort by neuro_deficit descending (True=1, False=0),
    # then by instability_score descending.
    # The `lambda` function creates a tuple for sorting: (True/False, score).
    # `reverse=True` sorts both keys in descending order.
    prioritized_list = sorted(patients, key=lambda p: (p['neuro_deficit'], p['instability_score']), reverse=True)

    print("Spine Surgery Prioritization:")
    print("----------------------------")
    print("The primary criterion for prioritization is the presence of a neurologic deficit.")
    print("The secondary criterion is the degree of mechanical instability.")
    print("\nFinal Triage Order (Highest to Lowest Priority):\n")

    final_order = []
    for i, patient in enumerate(prioritized_list):
        priority = i + 1
        print(f"Priority #{priority}: Patient {patient['id']}")
        print(f"  - Diagnosis: {patient['diagnosis']}")
        reason = "Active Neurologic Deficit (Surgical Emergency)" if patient['neuro_deficit'] else "High Degree of Mechanical Instability"
        if not patient['neuro_deficit'] and patient['instability_score'] < 2:
            reason = "Moderate Mechanical Instability"
        print(f"  - Reason for Priority: {reason}\n")
        final_order.append(str(patient['id']))
    
    print("The final priority order is Patient " + ", Patient ".join(final_order) + ".")
    print("This corresponds to the answer choice showing the patient order 3, 1, 2.")

prioritize_spine_patients()
<<<F>>>