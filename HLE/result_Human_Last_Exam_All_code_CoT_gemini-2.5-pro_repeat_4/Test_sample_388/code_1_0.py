def provide_counseling():
    """
    This function provides a pharmacy counseling recommendation based on the patient's medications.
    """
    
    # Patient's prescriptions
    atorvastatin_dose = 20
    junel_fe_dose = "1.5/30"
    fluoxetine_dose = 20

    # Counseling recommendation text
    recommendation = f"""
Hello Allison, I have your prescriptions ready here. I'd like to take a moment to discuss your medications.

I see you are picking up:
- Atorvastatin {atorvastatin_dose}mg
- Junel Fe {junel_fe_dose}mg
- Fluoxetine {fluoxetine_dose}mg

You mentioned that you took Excedrin for a headache, and I want to give you an important safety recommendation about that.

Your antidepressant, Fluoxetine {fluoxetine_dose}mg, can increase the risk of stomach bleeding. The aspirin in Excedrin also carries a similar risk. When you take them together, that risk of bleeding becomes much higher.

For your safety, I strongly recommend that you avoid using Excedrin or any other products containing aspirin or ibuprofen while you are taking Fluoxetine.

A much safer option for your headaches would be a product that only contains acetaminophen, like Tylenol.

Please let me know if you have any other questions.
"""

    print(recommendation)

provide_counseling()