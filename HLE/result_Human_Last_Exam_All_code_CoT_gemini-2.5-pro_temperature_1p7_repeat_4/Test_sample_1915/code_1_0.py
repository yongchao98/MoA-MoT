def solve_clinical_vignette():
    """
    This function analyzes the clinical case to determine the damaged anatomical structure.
    """

    # Step 1: List the patient's symptoms in the right eye.
    symptoms = [
        "No pupillary light reflex (controlled by CN III parasympathetics)",
        "Inability to adduct (medial rectus, innervated by CN III)",
        "Inability to depress (inferior rectus, innervated by CN III)",
        "Inability to elevate (superior rectus & inferior oblique, innervated by CN III)"
    ]

    # Step 2: Identify the common nerve responsible for these functions.
    conclusion_nerve = "All of these functions are controlled by Cranial Nerve III (Oculomotor Nerve)."

    # Step 3: Identify the anatomical location of Cranial Nerve III's nucleus.
    location_of_nucleus = "The nuclei for Cranial Nerve III (the oculomotor nucleus and the Edinger-Westphal nucleus) are located in the midbrain."

    # Step 4: Evaluate the options and conclude.
    reasoning = (
        "The patient's presentation is a classic right-sided Cranial Nerve III palsy.\n"
        "Given the answer choices, the damage must be to the anatomical structure that contains the origin of CN III.\n"
        "This structure is the midbrain."
    )

    final_answer = "E"

    # Print the reasoning and the final answer.
    print("Clinical Reasoning:")
    for symptom in symptoms:
        print(f"- {symptom}")
    print("\nConclusion about the Nerve:")
    print(conclusion_nerve)
    print("\nAnatomical Location:")
    print(location_of_nucleus)
    print("\nFinal Reasoning:")
    print(reasoning)
    print(f"\nTherefore, the correct answer is choice {final_answer}.")

solve_clinical_vignette()