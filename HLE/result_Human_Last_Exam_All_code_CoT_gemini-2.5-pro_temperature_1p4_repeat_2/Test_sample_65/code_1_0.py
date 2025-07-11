def prioritize_spine_patients():
    """
    Prioritizes spine surgery patients based on diagnosis and neurologic status.
    The function assigns a priority score where a lower score indicates higher urgency.
    """

    patients = [
        {
            "id": "Patient 1",
            "diagnosis": "Severe burst fracture of L2",
            "neuro_status": "No neurologic deficits",
            "priority": 0
        },
        {
            "id": "Patient 2",
            "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1",
            "neuro_status": "No neurologic deficits",
            "priority": 0
        },
        {
            "id": "Patient 3",
            "diagnosis": "Split fracture of L2",
            "neuro_status": "Mildly disordered pelvic functions",
            "priority": 0
        }
    ]

    # Assign priority scores based on clinical urgency
    # 1: Active neurologic compromise (surgical emergency)
    # 2: High-grade spinal instability (urgent)
    # 3: Lower-grade instability (less urgent)
    for p in patients:
        if "disordered pelvic functions" in p["neuro_status"]:
            p["priority"] = 1
            p["rationale"] = "Active neurologic compromise (cauda equina syndrome) is a surgical emergency."
        elif "Severe burst fracture" in p["diagnosis"]:
            p["priority"] = 2
            p["rationale"] = "High-grade instability poses a significant risk for future neurologic injury."
        else:
            p["priority"] = 3
            p["rationale"] = "Lower-grade instability without neurologic deficits presents the lowest immediate risk."

    # Sort patients by their priority score (ascending)
    sorted_patients = sorted(patients, key=lambda x: x["priority"])

    print("Surgical Prioritization from Highest to Lowest Priority:\n")
    for i, p in enumerate(sorted_patients, 1):
        print(f"Priority #{i}: {p['id']}")
        print(f"  Diagnosis: {p['diagnosis']}")
        print(f"  Neurologic Status: {p['neuro_status']}")
        print(f"  Rationale: {p['rationale']}\n")

    # Output the final ordered list
    final_order = ", ".join([p['id'] for p in sorted_patients])
    print(f"Final Prioritized Order: {final_order}")

prioritize_spine_patients()
<<<F>>>