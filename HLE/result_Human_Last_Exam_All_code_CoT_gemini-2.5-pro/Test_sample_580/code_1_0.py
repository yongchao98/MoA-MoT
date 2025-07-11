import textwrap

def diagnose_maneuver():
    """
    This function analyzes the clinical vignette to determine the correct
    physical exam maneuver.
    """

    # Step 1: Analyze the clinical information provided.
    patient_info = {
        "Complaint": "Right lower extremity pain in L4-S1 distribution (sciatica).",
        "Key Feature": "Pain is intensified by lying supine.",
        "Imaging": "X-ray imaging is unremarkable.",
        "Conclusion": "The symptoms point towards sciatica, but the pattern suggests a non-spinal cause like Piriformis Syndrome."
    }

    print("Step 1: Analyzing the Patient's Clinical Presentation")
    print("-----------------------------------------------------")
    for key, value in patient_info.items():
        print(f"- {key}: {value}")
    print("\n")

    # Step 2: Explain the rationale for Piriformis Syndrome.
    explanation = """
    Piriformis syndrome is a condition where the piriformis muscle, located deep in the buttock, compresses the sciatic nerve. This causes pain, tingling, and numbness along the path of the nerve, mimicking sciatica from a spinal cause. The piriformis muscle's primary function is to rotate the thigh bone outward (external rotation).
    """
    print("Step 2: Identifying the Most Likely Diagnosis to Test For")
    print("---------------------------------------------------------")
    print(textwrap.fill(explanation, width=80))
    print("\n")

    # Step 3: Analyze the physical exam setup.
    exam_setup = {
        "Patient Position": "Left decubitus (lying on the left side).",
        "Leg Tested": "Right leg (the symptomatic leg).",
        "Leg Position": "Extended.",
        "Test Type": "Resisted active motion (patient moves against physician's resistance)."
    }
    print("Step 3: Evaluating the Physical Exam Maneuver")
    print("----------------------------------------------")
    for key, value in exam_setup.items():
        print(f"- {key}: {value}")
    print("- Goal: To make the suspected muscle (piriformis) contract to see if it reproduces the pain.\n")


    # Step 4: Determine the correct action.
    actions = {
        "A. Abduction": "Primarily tests gluteus medius/minimus. Piriformis is a weak assistant.",
        "B. Adduction": "Tests the adductor muscles; relaxes the piriformis.",
        "C. Internal Rotation": "Opposite of piriformis function; would be a passive stretch, not a resisted contraction.",
        "D. External Rotation": "This is the PRIMARY function of the piriformis muscle. Resisting this action directly tests the muscle.",
        "E. Flexion": "Tests hip flexors (e.g., iliopsoas).",
        "F. Extension": "Tests gluteus maximus and hamstrings."
    }
    
    correct_action = "D"
    
    print("Step 4: Connecting the Action to the Diagnosis")
    print("-----------------------------------------------")
    print("To confirm piriformis syndrome via a resisted test, the patient must perform the muscle's primary action.")
    for action_key, desc in actions.items():
        if action_key.startswith(correct_action):
            print(f"-> {action_key}: {desc} <-- This is the correct maneuver.")
        else:
            print(f"- {action_key}: {desc}")
            
    print("\nFinal Conclusion:")
    print("Performing resisted external rotation of the extended hip will cause the piriformis to contract, which would compress the sciatic nerve and reproduce the patient's symptoms if piriformis syndrome is the cause.")
    print("The correct answer is D.")

diagnose_maneuver()