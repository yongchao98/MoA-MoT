def prioritize_patients():
    """
    Prioritizes spine surgery patients based on their diagnosis and neurological status.
    """
    patients = [
        {
            "id": 1,
            "name": "Patient 1",
            "diagnosis": "Severe burst fracture of L2, no neurologic deficits.",
            "priority": 0,
            "rationale": "Severe instability from a burst fracture poses a high risk of future neurological injury. This is a high priority for stabilization, but less urgent than an active neurological deficit."
        },
        {
            "id": 2,
            "name": "Patient 2",
            "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits.",
            "priority": 0,
            "rationale": "Instability from the fracture and mild slippage warrants surgery, but it is considered less mechanically unstable than a severe burst fracture, placing it at a lower priority."
        },
        {
            "id": 3,
            "name": "Patient 3",
            "diagnosis": "Split fracture of L2, with mildly disordered pelvic functions.",
            "priority": 0,
            "rationale": "The presence of 'disordered pelvic functions' indicates cauda equina syndrome, a neurologic emergency that requires immediate surgical decompression to prevent permanent damage. This is the top priority."
        }
    ]

    # Assign priority scores (higher is more urgent)
    for p in patients:
        if "disordered pelvic functions" in p["diagnosis"]:
            p["priority"] = 3  # Neurologic emergency
        elif "Severe burst fracture" in p["diagnosis"]:
            p["priority"] = 2  # High-grade instability
        elif "Compression fracture" in p["diagnosis"]:
            p["priority"] = 1  # Moderate instability

    # Sort patients by priority in descending order
    sorted_patients = sorted(patients, key=lambda x: x["priority"], reverse=True)

    print("Patient Prioritization from Highest to Lowest Priority:\n")
    for i, p in enumerate(sorted_patients):
        print(f"Priority #{i+1}: {p['name']} (ID: {p['id']})")
        print(f"  Diagnosis: {p['diagnosis']}")
        print(f"  Rationale: {p['rationale']}\n")
    
    final_order = [p['id'] for p in sorted_patients]
    print(f"The final prioritization order is: Patient {final_order[0]}, Patient {final_order[1]}, Patient {final_order[2]}")


prioritize_patients()
<<<F>>>