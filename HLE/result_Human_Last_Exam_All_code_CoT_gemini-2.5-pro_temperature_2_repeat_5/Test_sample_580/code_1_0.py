import sys

def solve_clinical_case():
    """
    This function analyzes the clinical vignette to determine the correct physical exam maneuver.
    """

    # Clinical Analysis:
    # 1. The patient's symptoms are pain in the L4-S1 dermatome, which is consistent with sciatica (sciatic nerve irritation).
    # 2. An unremarkable X-ray makes a bony cause (like a fracture or severe arthritis) less likely, pointing towards a soft-tissue or neuromuscular cause.
    # 3. Piriformis syndrome is a condition where the piriformis muscle compresses the sciatic nerve, causing sciatica. It is a strong candidate for the diagnosis.
    # 4. The piriformis muscle is a primary EXTERNAL rotator of the hip.
    # 5. Physical exam tests for piriformis syndrome aim to provoke symptoms by either stretching the muscle or having it contract against resistance.
    # 6. Stretching the piriformis muscle is achieved by performing the opposite action of its primary function.
    # 7. Therefore, to stretch the external rotator, the physician must perform INTERNAL ROTATION of the hip.
    # 8. This maneuver (passive forced internal rotation of the extended thigh) is a known clinical sign called Freiberg's sign, used to help diagnose piriformis syndrome.

    # Final Answer Determination
    answer_choice = "C"
    action_description = "Internal Rotation"

    # Print the logical deduction as the output
    print("Patient's symptoms suggest piriformis syndrome, where the piriformis muscle irritates the sciatic nerve.")
    print("The primary action of the piriformis muscle is external rotation of the hip.")
    print("To test for this condition, a physician can stretch the muscle to see if it reproduces the pain.")
    print("The action that stretches this external rotator muscle is internal rotation.")
    print(f"Therefore, the correct action to help confirm the diagnosis is: {answer_choice}. {action_description}")

# Execute the function to provide the solution.
solve_clinical_case()