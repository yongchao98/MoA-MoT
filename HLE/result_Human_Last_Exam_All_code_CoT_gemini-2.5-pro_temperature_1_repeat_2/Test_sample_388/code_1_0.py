def generate_counseling_recommendation():
    """
    Generates a pharmacist's counseling recommendation based on a patient's medications.
    """
    # Patient's medications from the scenario
    otc_med = "Excedrin"
    prescription1 = "Atorvastatin 20mg"
    prescription2 = "Junel Fe 1.5/30mg"
    prescription3 = "Fluoxetine 20mg"

    # The core counseling recommendation text
    recommendation = f"""
Hello Allison, I have your prescriptions ready here for Atorvastatin 20mg, Junel Fe 1.5/30mg, and Fluoxetine 20mg.

I noticed you mentioned taking Excedrin for your headache. I want to provide an important piece of counseling regarding that. Your new prescription, Fluoxetine 20mg, can increase the risk of stomach bleeding. The aspirin in most Excedrin products also carries a similar risk. When you take them together, that risk of bleeding becomes much higher.

For future headaches, I would recommend a safer alternative like acetaminophen (the active ingredient in Tylenol) instead of Excedrin.

It is a great choice for pain relief that doesn't have this specific interaction with your Fluoxetine 20mg. Please feel free to ask if you have any other questions about your medications.
"""
    print(recommendation)

generate_counseling_recommendation()