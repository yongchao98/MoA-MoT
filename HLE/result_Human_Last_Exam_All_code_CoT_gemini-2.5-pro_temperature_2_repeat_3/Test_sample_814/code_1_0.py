def solve_medical_case():
    """
    Analyzes a patient's clinical presentation and determines the best treatment option.
    """

    # Step 1: Analyze Patient's Profile
    symptoms = [
        "Chronic widespread pain (>1 year)",
        "Extreme fatigue",
        "Anxiety and depression",
        "Sleep issues",
        "Diminished cognitive ability ('fibro fog')",
        "Restless leg syndrome",
        "Paresthesia (numbness/tingling)"
    ]
    rule_outs = [
        "Normal thyroid function",
        "Negative for Rheumatoid Arthritis (RA)",
        "Negative for Lupus (SLE)",
        "Normal ESR (low inflammation)"
    ]

    print("--- Medical Case Analysis ---")
    print("\nStep 1: Patient Diagnosis")
    print("The patient's symptom profile, including chronic widespread pain, fatigue, sleep and mood disturbances, combined with the absence of inflammatory markers, strongly indicates a diagnosis of Fibromyalgia.")

    # Step 2: Evaluate Treatment Options
    print("\nStep 2: Evaluating Treatment Options for Fibromyalgia")

    print("\nAnalysis of Key Medications:")
    print(" - Duloxetine: An SNRI antidepressant, it is FDA-approved for fibromyalgia. It targets the core symptoms of pain as well as co-occurring anxiety and depression.")
    print(" - Gabapentin: An anticonvulsant, it is highly effective for neuropathic pain, which is a key component of fibromyalgia. It also effectively treats paresthesia, restless leg syndrome, and can improve sleep quality.")
    print(" - Ibuprofen (current med): An NSAID, which offers mild relief but does not address the central pain mechanism of fibromyalgia.")

    print("\nAnalysis of Combination Options:")
    print(" - A. Duloxetine+Gabapentin: This combination offers comprehensive, synergistic benefits. Duloxetine addresses pain and mood, while Gabapentin specifically targets the neuropathic pain, restless leg syndrome, and sleep issues. This covers the patient's full symptom profile.")
    print(" - B. Gabapentin alone: Helpful, but may not sufficiently address the patient's anxiety and depression.")
    print(" - C. Duloxetine alone: A strong first-line option, but may not be as effective as a combination for the patient's severe sleep issues and restless leg syndrome.")
    print(" - D. Cyclobenzaprine: A muscle relaxant mainly used as a sleep aid in fibromyalgia; not a comprehensive treatment.")
    print(" - E. Duloxetine+Acetaminophen: Adding acetaminophen offers little benefit over the patient's current ibuprofen use and doesn't address sleep or nerve symptoms.")
    print(" - F. Duloxetine+Cyclobenzaprine: A good combination for pain and sleep, but Gabapentin is superior to cyclobenzaprine for treating restless leg syndrome and paresthesia.")

    # Step 3: Conclusion
    final_choice = "A. Duloxetine+Gabapentin"
    reasoning = "This combination provides the most comprehensive treatment for the patient's multifaceted symptoms, addressing pain, mood (anxiety/depression), sleep disturbance, and specific neurological symptoms (restless leg syndrome, paresthesia) through two different mechanisms of action."

    print("\nStep 3: Final Conclusion")
    print(f"The best choice is {final_choice}.")
    print(f"Reasoning: {reasoning}")


if __name__ == "__main__":
    solve_medical_case()