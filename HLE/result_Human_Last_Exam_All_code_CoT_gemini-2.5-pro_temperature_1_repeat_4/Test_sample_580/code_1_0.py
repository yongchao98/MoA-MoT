def solve_clinical_case():
    """
    This script outlines the clinical reasoning to solve the patient case.
    It does not perform calculations but explains the logic step-by-step.
    """
    # Patient data from the prompt
    patient_age = 55
    pain_duration_months = 5
    pain_distribution = "L4-S1"
    bmi = 29

    # Step 1: Analyze the core problem.
    # The patient's symptoms (pain in the L4-S1 distribution) are characteristic of sciatica.
    # An unremarkable X-ray suggests a soft-tissue cause, not a bony one.
    print(f"Analysis of the {patient_age}-year-old patient with a {pain_duration_months}-month history of pain:")
    print(f"The L4-S1 pain pattern points towards sciatica. With a normal X-ray, we consider soft tissue causes like Piriformis Syndrome.")

    # Step 2: Understand the function of the piriformis muscle.
    # Piriformis Syndrome occurs when the piriformis muscle compresses the sciatic nerve.
    # The primary action of the piriformis muscle is to externally rotate the hip.
    print("\nThe key to the diagnosis is understanding the function of the muscle being tested.")
    print("The Piriformis muscle is a primary EXTERNAL ROTATOR of the hip.")

    # Step 3: Evaluate the physical exam maneuver.
    # The test described involves applying resistance to the extended right leg while the patient is lying on her left side.
    # This is a provocative test designed to reproduce symptoms by contracting a specific muscle.
    # Contracting the piriformis muscle against resistance would compress the sciatic nerve (if it's the cause) and reproduce the pain.
    print("\nThe physical exam is a provocative test. To confirm Piriformis Syndrome, the physician must make the muscle contract.")
    
    # Step 4: Conclude the correct action.
    # The action that causes the piriformis muscle to contract is external rotation.
    action = "D. External Rotation"
    print(f"\nConclusion: The correct action to provoke the piriformis and confirm the diagnosis is {action}.")

solve_clinical_case()