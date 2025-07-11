import sys

def solve_medical_case():
    """
    This function analyzes the clinical case and determines the correct physical exam maneuver.
    """
    
    # Step 1: Identify the likely diagnosis from the patient's symptoms.
    # The patient's symptoms (L4-S1 pain distribution without remarkable X-ray findings)
    # strongly suggest Piriformis Syndrome, where the piriformis muscle compresses the sciatic nerve.
    print("Clinical Analysis Steps:")
    print("1. The patient's presentation with radicular pain in the L4-S1 distribution and a normal X-ray points towards a soft tissue cause of sciatica, such as Piriformis Syndrome.")
    
    # Step 2: Understand the function of the piriformis muscle.
    # The piriformis muscle's primary action is the external rotation of the femur at the hip joint.
    print("2. The main function of the piriformis muscle, located in the buttock, is to externally rotate the hip.")
    
    # Step 3: Determine the correct provocative maneuver.
    # To confirm Piriformis Syndrome, a physician performs a provocative test.
    # Having the patient activate the muscle against resistance will reproduce the pain if the muscle is inflamed or spasming.
    # The test is performed with the patient in the left decubitus position (lying on the non-symptomatic side).
    # Resisting the primary action of the muscle is the most direct test.
    print("3. To confirm Piriformis Syndrome, a physician will apply resistance while the patient actively engages the muscle. Given its function, this action is External Rotation of the hip.")
    
    # Step 4: Conclude with the correct answer choice.
    print("4. Therefore, asking the patient to externally rotate their extended right leg against resistance will provoke pain and help confirm the diagnosis of Piriformis Syndrome.")

    final_answer = "D"
    sys.stdout.write(f"\nFinal Answer:\n<<<{final_answer}>>>\n")

solve_medical_case()