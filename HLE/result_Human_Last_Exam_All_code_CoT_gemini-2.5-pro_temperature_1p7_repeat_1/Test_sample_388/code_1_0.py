def pharmacy_counseling():
    """
    Generates a counseling recommendation for Allison based on her medications and symptoms.
    """
    # Medications and symptoms
    prescriptions = {
        "Atorvastatin": "20mg",
        "Junel Fe": "1.5/30mg",
        "Fluoxetine": "20mg"
    }
    otc_medication = "Excedrin (contains Aspirin)"
    symptoms = "New headache with an unusual smell (potential aura)"

    # Formulate counseling points
    # Primary Concern: New headache while on a combined oral contraceptive (Junel Fe)
    # This is a red flag for a potential thromboembolic event (e.g., stroke).
    primary_recommendation = (
        f"The most important concern is your new headache. "
        f"When taking a birth control pill like Junel Fe {prescriptions['Junel Fe']}, "
        f"a new or worsening headache, especially with an unusual symptom like a strange smell, "
        f"can be a warning sign of a serious side effect like a blood clot. "
        f"You must contact your doctor immediately to discuss this."
    )

    # Secondary Concern: Interaction between Fluoxetine (SSRI) and Aspirin (NSAID in Excedrin)
    secondary_recommendation = (
        f"Also, please be aware that taking your antidepressant, Fluoxetine {prescriptions['Fluoxetine']}, "
        f"with aspirin, which is in the Excedrin you took, increases the risk of bleeding. "
        f"Once you have spoken to your doctor about the headache, you should discuss safer pain relief options."
    )

    # Final recommendation printout
    final_message = (
        "Here is the counseling recommendation based on your information:\n\n"
        "Pharmacist to Allison:\n\"" +
        primary_recommendation +
        "\n\n" +
        secondary_recommendation +
        "\""
    )

    print(final_message)

pharmacy_counseling()