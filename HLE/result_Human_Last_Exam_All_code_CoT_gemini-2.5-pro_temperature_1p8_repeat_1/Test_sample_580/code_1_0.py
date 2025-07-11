import textwrap

def solve_clinical_case():
    """
    This function analyzes a clinical vignette and determines the correct diagnostic maneuver.
    The reasoning is printed step-by-step.
    """

    print("### Clinical Case Analysis ###\n")

    # Step 1: Analyze the patient's presentation
    print("Step 1: Analyzing the Patient's Presentation")
    print("-" * 40)
    analysis_text = """
    The patient presents with sciatica-like symptoms (pain in the L4-S1 dermatome) with a normal X-ray. This suggests a soft-tissue cause for nerve compression, such as a herniated disc or, importantly, Piriformis Syndrome. Piriformis Syndrome occurs when the piriformis muscle irritates or compresses the sciatic nerve.
    """
    print(textwrap.fill(analysis_text, width=80))
    print("\n")

    # Step 2: Understand the physical exam for Piriformis Syndrome
    print("Step 2: Understanding the Physical Exam")
    print("-" * 40)
    exam_text = """
    The physician is examining the patient in the left decubitus (side-lying) position to test the right hip. To diagnose Piriformis Syndrome, a clinician will perform a maneuver to stretch the piriformis muscle, which in turn compresses the sciatic nerve and should reproduce the patient's specific pain if the syndrome is present.
    """
    print(textwrap.fill(exam_text, width=80))
    print("\n")
    
    # Step 3: Determine the correct maneuver
    print("Step 3: Determining the Correct Maneuver")
    print("-" * 40)
    maneuver_text = """
    The primary action of the piriformis muscle is EXTERNAL rotation of the hip. Therefore, to stretch the muscle, the physician must perform the opposite motion. The opposite motion is INTERNAL rotation. By passively forcing the hip into internal rotation, the piriformis muscle is stretched, and if it is the source of the sciatic nerve compression, the patient's pain will be reproduced. This is a key component of tests like the Freiberg test or the FAIR (Flexion, Adduction, Internal Rotation) test.
    """
    print(textwrap.fill(maneuver_text, width=80))
    print("\n")

    # Step 4: Final Conclusion
    print("Step 4: Conclusion")
    print("-" * 40)
    print("The action that will stretch the piriformis muscle and confirm the diagnosis is Internal Rotation.")
    print("\n")
    
    final_answer = "C"
    print(f"The final answer is {final_answer}.")


# Execute the analysis
solve_clinical_case()

# Final answer in the required format
print("<<<C>>>")