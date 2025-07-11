def generate_pharmacist_recommendation():
    """
    This function formulates and prints a pharmacist's counseling recommendation
    based on a patient's medication profile and symptoms.
    """

    # Patient's prescriptions with dosages from the prompt
    atorvastatin = "Atorvastatin 20mg"
    junel_fe = "Junel Fe 1.5/30mg"
    fluoxetine = "Fluoxetine 20mg"

    # OTC medication used by the patient
    otc_medication = "Excedrin"

    # Print the counseling points in a clear, user-friendly format
    print("Based on the information provided, here is a potential counseling recommendation from the pharmacist:\n")
    print("-------------------------------------------------------------------------------------------------")
    print(f"Hello Allison, I have your prescriptions for {atorvastatin}, {junel_fe}, and {fluoxetine} ready.")
    print(f"I understand you took some {otc_medication} for a headache you experienced recently. I'd like to give you some important information about your medicines.\n")
    
    # Counseling point 1: The serious side effect of the oral contraceptive
    print("1. Your birth control, Junel Fe 1.5/30mg, has some rare but serious side effects.")
    print("   A new or severe headache can be a warning sign. It is very important that you follow up with your doctor about this new headache to make sure it is not a sign of a more serious issue.\n")
    
    # Counseling point 2: The drug interaction
    print("2. Your antidepressant, Fluoxetine 20mg, can interact with the aspirin found in Excedrin.")
    print("   Taking them together increases the risk of stomach or intestinal bleeding. For managing simple headaches in the future, it would be safer for you to use a product containing only acetaminophen, like Tylenol.\n")

    print("So, the key takeaways are to please contact your doctor about the headache and to use acetaminophen instead of Excedrin for pain relief moving forward.")
    print("-------------------------------------------------------------------------------------------------")

generate_pharmacist_recommendation()