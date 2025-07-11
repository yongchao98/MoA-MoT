def provide_counseling():
    """
    This function provides a pharmacist's counseling recommendation based on the patient's medications and symptoms.
    """
    
    # Medications
    otc_med = "Excedrin (contains Aspirin)"
    prescription1 = "Fluoxetine 20mg"
    prescription2 = "Junel Fe 1.5/30mg"
    prescription3 = "Atorvastatin 20mg"
    
    # Main counseling points
    interaction_med1 = "Fluoxetine"
    interaction_med2 = "Aspirin (in Excedrin)"
    headache_risk_med = "Junel Fe (your oral contraceptive)"
    
    print("Pharmacist's Counseling Recommendation for Allison:")
    print("--------------------------------------------------")
    print(f"Hello, I have your prescriptions ready. I noticed you picked up some {otc_med} for a headache.")
    print("\nI want to make you aware of two important points regarding your medicines:\n")
    
    # Point 1: Drug Interaction
    print(f"1. Increased Bleeding Risk: Your prescription for {interaction_med1} can increase the risk of bleeding.")
    print(f"   The {interaction_med2} also thins the blood. Taking these two medications together significantly increases the risk of stomach bleeding.")
    print("   For future headaches, it would be safer to use a product with just acetaminophen (like Tylenol) instead.\n")

    # Point 2: Symptom Warning
    print(f"2. Follow Up on Your Headache: It is very important that you contact your doctor about this new headache.")
    print(f"   New or severe headaches can be a rare but serious warning sign for patients taking {headache_risk_med}.")
    print("   Your doctor needs to rule out anything serious, so please give their office a call soon to let them know.\n")

    print("Please let me know if you have any other questions.")

# Execute the function to print the counseling recommendation
provide_counseling()