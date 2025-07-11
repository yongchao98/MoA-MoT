def solve_clinical_case():
    """
    This script analyzes the clinical vignette and determines the correct physical exam maneuver.
    """
    
    # Step 1: Analyze the patient's key symptoms and signs.
    # The patient has right-sided L4-S1 radicular pain with a normal X-ray.
    # This suggests a soft-tissue cause of sciatic nerve compression, such as Piriformis Syndrome.
    print("Patient Presentation Analysis:")
    print("The patient's symptoms point towards sciatica without a clear bony cause, making Piriformis Syndrome a likely diagnosis.")
    print("-" * 20)

    # Step 2: Understand the function of the piriformis muscle.
    # The piriformis muscle's primary action is external rotation of the femur at the hip joint.
    print("Anatomical Consideration:")
    print("The piriformis muscle is a primary external rotator of the hip. Piriformis syndrome occurs when this muscle compresses the sciatic nerve.")
    print("-" * 20)

    # Step 3: Analyze the described physical exam maneuver.
    # The physician is applying resistance to the patient's leg movement. This is a resisted muscle test.
    # The goal is to make the suspected muscle (piriformis) contract to see if it reproduces the pain.
    print("Physical Exam Maneuver Analysis:")
    print("The test described is a resisted muscle test designed to provoke symptoms by contracting a specific muscle.")
    print("To test the piriformis muscle, the patient must be asked to perform its primary action against resistance.")
    print("-" * 20)

    # Step 4: Conclude the correct action.
    # Contracting the piriformis muscle is achieved through external rotation.
    # Therefore, resisted external rotation would reproduce the pain and help confirm the diagnosis.
    final_action = "External Rotation"
    answer_choice = "D"
    
    print("Conclusion:")
    print(f"The action that contracts the piriformis muscle is {final_action}.")
    print(f"Performing resisted {final_action} will reproduce the patient's pain if they have Piriformis Syndrome.")
    print(f"This corresponds to answer choice {answer_choice}.")

solve_clinical_case()