def triage_spine_patients():
    """
    Prioritizes three spine surgery patients based on their clinical presentation.
    """
    # Patient data with assigned scores for neurological status and instability
    # Scoring Rationale:
    # Neurological Score: Disordered pelvic function (Cauda Equina sign) is a surgical emergency (10 points). No deficit is 0.
    # Instability Score: Split fracture (highly unstable) = 8. Severe burst fracture (highly unstable) = 7. Compression + mild slip (less unstable) = 5.
    patients = [
        {
            "id": "Patient 1",
            "diagnosis": "Severe burst fracture of L2, no neurologic deficits",
            "neuro_score": 0,
            "instability_score": 7
        },
        {
            "id": "Patient 2",
            "diagnosis": "Compression fracture of L2 with mild traumatic spondylolisthesis of L1, no neurologic deficits",
            "neuro_score": 0,
            "instability_score": 5
        },
        {
            "id": "Patient 3",
            "diagnosis": "Split fracture of L2, with mildly disordered pelvic functions",
            "neuro_score": 10,
            "instability_score": 8
        }
    ]

    # Calculate total priority score for each patient
    for p in patients:
        p["priority_score"] = p["neuro_score"] + p["instability_score"]

    # Sort patients in descending order of priority score
    sorted_patients = sorted(patients, key=lambda x: x["priority_score"], reverse=True)

    print("Spine Surgery Triage Prioritization:\n")
    print("Priority is determined by neurological status and spinal stability.")
    print("A neurological deficit, especially involving pelvic function, is a surgical emergency.\n")

    for i, p in enumerate(sorted_patients):
        print(f"Priority #{i+1}: {p['id']}")
        print(f"  Diagnosis: {p['diagnosis']}")
        print(f"  Reasoning: Priority score of {p['priority_score']} (Neuro: {p['neuro_score']}, Instability: {p['instability_score']})\n")

    final_order = [p['id'].replace('Patient ', '') for p in sorted_patients]
    print(f"Final Prioritization (Highest to Lowest): Patient {final_order[0]}, Patient {final_order[1]}, Patient {final_order[2]}")

if __name__ == '__main__':
    triage_spine_patients()
<<<F>>>