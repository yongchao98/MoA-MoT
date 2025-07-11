def provide_counseling_recommendation():
    """
    Analyzes the patient's medications and provides a key counseling point.
    """
    
    # Patient's prescriptions
    medication_1 = "Atorvastatin 20mg"
    medication_2 = "Junel Fe 1.5/30mg"
    medication_3 = "Fluoxetine 20mg"
    
    # Over-the-counter medication taken
    otc_med = "Excedrin"
    
    # Explanation of the interaction
    explanation = (
        f"The most important counseling point is the potential interaction between her new prescription, {medication_3}, "
        f"and the {otc_med} she took for her headache.\n\n"
        f"Fluoxetine is a Selective Serotonin Reuptake Inhibitor (SSRI). {otc_med} contains Aspirin, which is a "
        "Nonsteroidal Anti-Inflammatory Drug (NSAID).\n\n"
        "When an SSRI like Fluoxetine is taken with an NSAID like Aspirin, there is a significantly "
        "increased risk of bleeding, particularly stomach and intestinal bleeding."
    )
    
    # Counseling recommendation
    recommendation = (
        "Recommendation: You should advise Allison about this increased bleeding risk. "
        "Recommend that she avoid using Excedrin or other NSAIDs (such as ibuprofen or naproxen) "
        "while taking Fluoxetine. For future headaches, suggest a safer alternative like acetaminophen (Tylenol), "
        "which does not carry the same risk of interaction."
    )
    
    # Print the full counseling point
    print("--- Pharmacy Counseling Recommendation ---")
    print(explanation)
    print("\n" + recommendation)

provide_counseling_recommendation()