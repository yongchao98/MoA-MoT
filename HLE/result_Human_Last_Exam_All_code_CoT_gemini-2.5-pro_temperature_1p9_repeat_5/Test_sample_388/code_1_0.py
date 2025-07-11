def provide_counseling():
    """
    This function prints a pharmacist's counseling recommendation based on
    the patient's medications.
    """
    patient_name = "Allison"
    otc_drug = "Excedrin"
    prescription_drug = "Fluoxetine"
    otc_ingredient = "aspirin"
    safer_alternative = "acetaminophen (Tylenol)"
    interaction_risk = "increase the risk of bleeding, especially in the stomach"

    print("--- Pharmacist Counseling Recommendation ---")
    print(f"\nHello {patient_name}, I have your prescriptions ready.")
    print(f"I noticed you are purchasing {otc_drug} for your headache.")
    print(f"\nI need to give you an important piece of advice regarding your medications.")
    print(f"One of your prescriptions, {prescription_drug}, can interact with the {otc_ingredient} found in {otc_drug}.")
    print(f"Taking these two medications together can {interaction_risk}.")
    print(f"\nFor your headache, a safer option would be a product containing just {safer_alternative}, as it doesn't have this interaction.")
    print("\nPlease consider using the safer alternative while you are taking Fluoxetine.")
    print("\n--- End of Recommendation ---")

provide_counseling()