def solve_medical_puzzle():
    """
    This script analyzes a medical case to identify the correct diagnostic maneuver.
    The logic focuses on the differential diagnosis of sciatica and the specific function
    of the piriformis muscle.
    """

    print("### Clinical Analysis Step-by-Step ###")
    print("1. Patient's symptoms are consistent with sciatica (pain in the L4-S1 nerve root distribution). Given the unremarkable X-ray, the cause is likely soft tissue-related, such as a herniated disc or, importantly, piriformis syndrome.")
    print("2. Piriformis syndrome occurs when the piriformis muscle compresses or irritates the sciatic nerve. Therefore, a key diagnostic test involves provoking the piriformis muscle.")
    print("3. The primary anatomical function of the piriformis muscle is to externally rotate the hip when the leg is extended.")
    print("4. The patient is in the left decubitus position (lying on their left side), which places the symptomatic right leg on top. The physician applies resistance to the extended right leg, which indicates a test of active muscle contraction.")
    print("5. To test the piriformis muscle specifically, the patient must be asked to perform its primary function against resistance. In this position, that action is external rotation of the right leg.")
    print("6. If performing external rotation against resistance reproduces the patient's characteristic L4-S1 pain, it strongly supports the diagnosis of piriformis syndrome over other causes like a disc herniation.")

    # The prompt requires an "equation" with numbers from the case.
    # We will format the diagnostic logic as an equation to meet this requirement.
    print("\n### Final Diagnostic Equation ###")
    print("For a 55-year-old patient with a 5-month history of pain in the L4-S1 distribution:")
    print("1 maneuver (External Rotation) + 1 source of resistance -> 1 confirmation of Piriformis Syndrome")


solve_medical_puzzle()