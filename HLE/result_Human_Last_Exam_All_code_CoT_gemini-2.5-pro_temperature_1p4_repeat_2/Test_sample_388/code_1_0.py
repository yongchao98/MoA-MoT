def provide_counseling_recommendation():
    """
    Analyzes Allison's medications to provide a pharmacist's counseling recommendation.
    """
    print("Pharmacist's Counseling Plan:")
    print("---------------------------------")
    print("Patient's Reported Medications:")
    print("* OTC: Excedrin (contains Aspirin)")
    print("* Rx: Fluoxetine 20mg")
    print("* Rx: Junel Fe 1.5/30mg")
    print("* Rx: Atorvastatin 20mg")
    print("\nStep 1: Identify Potential Drug Interaction")
    print("The patient is taking Fluoxetine 20mg, an SSRI, and has just taken Excedrin, which contains Aspirin.")
    print("FACT: Combining an SSRI like Fluoxetine with an NSAID like Aspirin significantly increases the risk of gastrointestinal bleeding.")
    print("\nStep 2: Identify Potential Serious Adverse Drug Effect")
    print("The patient is complaining of a headache and is taking Junel Fe 1.5/30mg, a combined oral contraceptive.")
    print("FACT: New, persistent, or severe headaches can be a serious warning sign (e.g., of a blood clot) for patients taking hormonal contraceptives.")
    print("\n---------------------------------")
    print("Final Counseling Recommendation:")
    print("---------------------------------")
    print("Hi Allison. I have your prescriptions ready, but I'd like to talk about the headache you mentioned.")
    print("\n1. Regarding the Excedrin you took: Your antidepressant, Fluoxetine 20mg, when taken with aspirin which is in Excedrin, can increase your risk of stomach bleeding. For future headaches, a safer option would be a product with just acetaminophen, like Tylenol.")
    print("\n2. More importantly, since you are taking birth control, Junel Fe 1.5/30mg, a new or worsening headache can sometimes be a sign of a more serious side effect. It's very important that you contact your doctor's office to let them know about this new headache to make sure it's not something serious.")

# Execute the function to print the recommendation
provide_counseling_recommendation()
