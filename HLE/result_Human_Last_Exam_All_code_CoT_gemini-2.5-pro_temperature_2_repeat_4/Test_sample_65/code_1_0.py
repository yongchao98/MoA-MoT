def triage_patients():
    """
    Prioritizes three spine surgery patients based on their clinical presentation.

    The prioritization is based on established neurosurgical principles:
    1. Active neurologic compromise (especially cauda equina syndrome) is a surgical emergency.
    2. Severe mechanical instability without deficits is highly urgent to prevent future deficits.
    3. Less severe instability without deficits is the lowest priority among these urgent cases.
    """
    
    # Patient data with assigned priority and reason
    # Priority 1 is the highest urgency.
    patients = [
        {
            "id": 1,
            "condition": "Severe burst fracture of L2, no neurologic deficits.",
            "priority": 2,
            "reason": "High priority due to severe mechanical instability and risk of future neurologic compromise."
        },
        {
            "id": 2,
            "condition": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.",
            "priority": 3,
            "reason": "Lowest priority of the three due to less severe instability (mild listhesis) and no neurologic deficits."
        },
        {
            "id": 3,
            "condition": "Split fracture of L2, with mildly disordered pelvic functions.",
            "priority": 1,
            "reason": "Top priority due to signs of cauda equina syndrome (disordered pelvic functions), a neurosurgical emergency."
        }
    ]

    # Sort patients based on the assigned priority number
    sorted_patients = sorted(patients, key=lambda x: x['priority'])

    print("Surgical Priority List:")
    print("=======================")
    for i, patient in enumerate(sorted_patients, 1):
        print(f"Priority #{i}: Patient {patient['id']}")
        print(f"  Condition: {patient['condition']}")
        print(f"  Reasoning: {patient['reason']}\n")
    
    # The final prioritized order is Patient 3, then Patient 1, then Patient 2.
    final_order = ", ".join([f"Patient {p['id']}" for p in sorted_patients])
    print(f"Final Prioritization Order: {final_order}")
    print("\nThis corresponds to Answer Choice F.")
    print("The final order is based on Patient 3 > Patient 1 > Patient 2")


triage_patients()
<<<F>>>