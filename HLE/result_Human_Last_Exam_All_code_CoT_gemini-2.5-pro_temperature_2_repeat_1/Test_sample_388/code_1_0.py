def provide_counseling():
    """
    Analyzes Allison's medications to provide a key counseling recommendation.
    """
    # Patient's prescriptions and recent OTC medication
    prescription_1 = "Atorvastatin 20mg"
    prescription_2 = "Junel Fe 1.5/30mg"
    prescription_3 = "Fluoxetine 20mg"
    otc_medication = "Excedrin (contains Aspirin)"

    # Explanation of the interaction
    print("Pharmacist's Analysis:")
    print("---------------------------------")
    print(f"Patient is taking {prescription_3}, which is an SSRI (Selective Serotonin Reuptake Inhibitor).")
    print(f"Patient also recently took {otc_medication}, which contains an NSAID (Nonsteroidal Anti-Inflammatory Drug).")
    print("\nInteraction Alert:")
    print("Both SSRIs and NSAIDs can increase the risk of bleeding, especially in the stomach.")
    print("Taking them together significantly increases this risk.")
    
    # Formulate the final recommendation
    print("\nCounseling Recommendation:")
    print("---------------------------------")
    recommendation = (
        "Because you are taking Fluoxetine, you should avoid using NSAIDs like the aspirin in Excedrin "
        "due to an increased risk of stomach bleeding. For headaches, a safer choice would be a medication with "
        "just acetaminophen, like Tylenol. Please speak with your doctor if your headaches continue."
    )
    print(recommendation)

provide_counseling()