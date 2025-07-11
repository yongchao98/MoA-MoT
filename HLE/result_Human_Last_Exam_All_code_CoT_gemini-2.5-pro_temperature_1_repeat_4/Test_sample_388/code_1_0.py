def pharmacist_counseling():
    """
    This function provides a pharmacist's counseling recommendation based on a patient's medications.
    """
    
    # Patient's prescriptions
    atorvastatin_dose = 20
    junel_fe_dose = "1.5/30"
    fluoxetine_dose = 20

    # OTC medication taken
    otc_medication = "Excedrin (contains Aspirin)"

    # The counseling recommendation
    print("Pharmacist Counseling Recommendation for Allison:")
    print("-" * 45)
    
    print(f"I see you are picking up your Atorvastatin {atorvastatin_dose}mg, Junel Fe, and Fluoxetine {fluoxetine_dose}mg today.")
    print(f"You mentioned taking Excedrin for a headache, and I want to share some important information.")
    
    print("\nPrimary Recommendation (Safety Interaction):")
    print(f"Your antidepressant, Fluoxetine, can interact with the aspirin found in Excedrin.")
    print("Taking these two medications together significantly increases the risk of stomach or intestinal bleeding.")
    print("For future headaches, a safer alternative would be a product containing only acetaminophen, such as Tylenol.")

    print("\nSecondary Recommendation (Monitoring):")
    print("With your birth control, Junel Fe, it is always a good idea to let your doctor know if you start having new or worsening headaches.")

    print("\nPlease let me know if you have any questions.")
    print("-" * 45)

pharmacist_counseling()
<<<Avoid taking Excedrin (which contains aspirin) with Fluoxetine due to an increased risk of stomach bleeding. A safer choice for headaches would be a product containing only acetaminophen, like Tylenol.>>>