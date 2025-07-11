def evaluate_trauma_patient():
    """
    This function analyzes the patient's clinical data to determine the
    first-line treatment for hemorrhagic shock.
    """
    # Patient Vitals and Key Findings
    heart_rate = 160  # beats per minute
    blood_pressure_systolic = 60 # mmHg (assuming 60/40)
    blood_pressure_diastolic = 40 # mmHg
    hemoglobin = 6  # gm/dL
    is_bleeding_profusely = True
    mental_status = "disoriented"
    skin = "cold and clammy"

    # Clinical Analysis
    is_tachycardic = heart_rate > 100
    is_hypotensive = blood_pressure_systolic < 90
    has_severe_anemia = hemoglobin < 7
    has_poor_perfusion = (mental_status == "disoriented" and skin == "cold and clammy")

    # Diagnosis based on findings
    diagnosis = None
    if is_tachycardic and is_hypotensive and has_severe_anemia and is_bleeding_profusely and has_poor_perfusion:
        diagnosis = "Hemorrhagic Shock"

    print(f"Patient Diagnosis: {diagnosis}")
    print("-------------------------------------")
    print("Rationale for Treatment:")
    print("The patient shows clear signs of hemorrhagic shock due to massive blood loss.")
    print("Immediate priority is to restore circulating volume to stabilize blood pressure and perfuse vital organs.")
    print("This is achieved with rapid intravenous infusion of isotonic crystalloid fluids.")
    print("-------------------------------------")

    # Treatment Options
    options = {
        'A': "Lay down the person and elevate legs along with CPR",
        'B': "Administer anticlotting medicine such as aspirin or heparin",
        'C': "Intravenous resuscitation of normal saline or Ringer's lactate",
        'D': "Intravenous resuscitation of normal saline",
        'E': "Intravenous resuscitation of normal saline with fructose"
    }

    # Selecting the best treatment
    correct_choice = None
    reasoning = ""
    if diagnosis == "Hemorrhagic Shock":
        correct_choice = 'C'
        reasoning = (
            "Option C is the most appropriate first-line treatment.\n"
            "Normal Saline and Ringer's Lactate are the standard isotonic crystalloids for initial volume resuscitation in trauma.\n"
            "Option A is incorrect because CPR is not indicated (patient has a pulse).\n"
            "Option B is dangerous as it would worsen bleeding.\n"
            "Options D and E are less ideal; D is less comprehensive than C, and E (sugar solution) is not used for initial shock resuscitation."
        )

    print(f"Recommended Treatment Option: {correct_choice}")
    print(f"Justification: {options[correct_choice]}")
    print("\nDetailed Reasoning:")
    print(reasoning)

    # Final Answer Block
    print("\n<<<C>>>")


if __name__ == "__main__":
    evaluate_trauma_patient()