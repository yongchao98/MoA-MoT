def provide_counseling():
    """
    Analyzes the patient's medications and provides a key counseling point.
    """

    # Patient's prescriptions
    prescriptions = {
        "Atorvastatin": "20mg",
        "Junel Fe": "1.5/30mg",
        "Fluoxetine": "20mg"
    }

    # Patient's OTC medication
    otc_medication = "Excedrin (contains Aspirin)"

    # The key interaction
    interacting_prescription = "Fluoxetine"
    interacting_otc_ingredient = "Aspirin"
    risk = "an increased risk of stomach and gastrointestinal bleeding"

    # The safer alternative
    alternative = "Acetaminophen (Tylenol)"

    recommendation = f"""
Hello Allison, I have your prescriptions ready. I overheard you mention that you took Excedrin for a headache and that you are also picking up your prescription for {interacting_prescription} {prescriptions[interacting_prescription]}.

It's very important for you to know that combining {interacting_prescription} with {interacting_otc_ingredient}, which is an ingredient in Excedrin, can cause {risk}.

For future headaches, I would recommend a safer alternative that doesn't have this interaction, such as plain {alternative}.

You should also watch for any signs of stomach bleeding, like severe stomach pain or dark, tarry stools, and contact your doctor if you experience them. If your headaches continue, it would be a good idea to discuss them with your doctor.
"""

    print(recommendation)

provide_counseling()