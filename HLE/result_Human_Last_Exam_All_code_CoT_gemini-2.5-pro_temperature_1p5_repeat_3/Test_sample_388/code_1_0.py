def provide_counseling():
    """
    This function generates a pharmacist's counseling recommendation
    based on Allison's medications and symptoms.
    """

    # Prescribed medications with dosages
    atorvastatin_dose = 20
    junel_fe_dose = "1.5/30"
    fluoxetine_dose = 20

    # The counseling recommendation string
    recommendation = f"""
Hello Allison, I have your prescriptions for Atorvastatin {atorvastatin_dose}mg, Junel Fe {junel_fe_dose}mg, and Fluoxetine {fluoxetine_dose}mg ready for you.

I noticed you were buying Excedrin for your headache, and I want to share some important information. The Fluoxetine {fluoxetine_dose}mg you are taking can interact with the aspirin component in Excedrin. When taken together, they can increase the risk of stomach upset and bleeding. For your headache, a safer alternative would be a product containing only acetaminophen, such as Tylenol.

Additionally, new or worsening headaches can sometimes be a side effect of oral contraceptives like your Junel Fe {junel_fe_dose}mg. It would be a good idea to mention these new headaches to your doctor, especially since you experienced an unusual smell beforehand. Your doctor will want to rule out any connection to your medication.
"""

    print(recommendation)

provide_counseling()