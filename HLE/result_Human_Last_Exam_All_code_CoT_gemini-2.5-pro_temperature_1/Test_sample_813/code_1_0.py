import textwrap

def solve_clinical_case():
    """
    This function analyzes the clinical case step-by-step to determine the root cause.
    """
    # Step 1: Analyze the patient's initial symptoms and history.
    initial_symptoms = ["agitation", "difficulty falling asleep", "increase in sexual activities (hypersexuality)", "careless spending"]
    family_history = "mood disorders"
    conclusion1 = "The initial symptoms and family history strongly suggest a manic episode, likely due to Bipolar Disorder."

    # Step 2: Identify the likely medication prescribed.
    treatment_for_mania = "A mood stabilizer, with Lithium being a very common and first-line choice."
    conclusion2 = "The patient was most likely prescribed Lithium to treat the manic episode."

    # Step 3: Analyze the new symptom that appeared after the prescription.
    new_symptom = "decreased interest in having sex (sexual dysfunction)"
    timeline = "The new symptom started *after* the medication was prescribed."
    conclusion3 = "This timeline points towards a medication side effect as the cause of the sexual dysfunction."

    # Step 4: Connect the likely medication (Lithium) to the new symptom.
    lithium_side_effect = "Lithium is well-known to cause hypothyroidism."
    hypothyroidism_symptom = "Hypothyroidism is a common cause of decreased libido and sexual dysfunction."
    conclusion4 = "Therefore, the most logical chain of events is: Bipolar Disorder -> Lithium prescription -> Lithium-induced hypothyroidism -> Sexual dysfunction."

    # Step 5: Evaluate the other answer choices based on the timeline.
    distractor_info = "30-year work history of metal smelting."
    distractor_analysis = """The heavy metal options (Lead, Arsenic, Mercury, Manganese) are less likely. While they can cause health issues, they do not explain the initial manic episode or the specific timing of the sexual dysfunction, which appeared only after the new medication was started."""

    # Final Conclusion
    final_conclusion = "The underlying root cause that connects the treatment for the initial condition to the final symptom is Lithium-induced hypothyroidism."
    answer_choice = "A"

    # Print the reasoning
    print("Clinical Reasoning Steps:")
    print("========================")
    print(f"1. Initial Diagnosis: {conclusion1}")
    print(f"2. Likely Prescription: {conclusion2}")
    print(f"3. Symptom Timeline: {conclusion3}")
    print(f"4. Connecting Treatment to Symptom: {conclusion4}")
    print(f"5. Evaluating Other Options: {textwrap.fill(distractor_analysis, 80)}")
    print("\nFinal Conclusion:")
    print("-----------------")
    print(f"{final_conclusion} This makes choice 'A' the best fit.")
    print("\nFinal Answer Choice:")
    print(f"'{answer_choice}'")

solve_clinical_case()
<<<A>>>