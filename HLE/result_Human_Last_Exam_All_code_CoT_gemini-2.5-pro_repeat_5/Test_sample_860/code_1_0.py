def assess_patient_and_recommend_treatment(hr, sbp, dbp, hb):
    """
    Assesses patient's condition based on vitals and recommends first-line treatment.
    """
    print("Patient's Key Clinical Data:")
    print(f"Heart Rate: {hr} bpm")
    print(f"Blood Pressure: {sbp}/{dbp} mmHg")
    print(f"Hemoglobin: {hb} gm/dL")
    print("-" * 30)

    # Check for signs of severe hemorrhagic shock
    if hr > 120 and sbp < 90 and hb < 7:
        print("Analysis: The patient's vital signs indicate severe hemorrhagic shock.")
        print("The combination of:")
        print(f" - Severe Tachycardia (HR of {hr})")
        print(f" - Profound Hypotension (BP of {sbp}/{dbp})")
        print(f" - Critical Anemia (Hemoglobin of {hb})")
        print("points to massive blood volume loss requiring immediate intervention.")
        print("\nPrimary Goal: Rapidly restore circulating volume to stabilize the patient.")
        print("\nRecommended First-Line Treatment:")
        print("Initiate aggressive intravenous (IV) fluid resuscitation with an isotonic crystalloid solution.")
        print("This makes option C the most appropriate choice.")

    else:
        print("The patient's condition requires further assessment.")

# Patient data from the case
heart_rate = 160
systolic_bp = 40
diastolic_bp = 60
hemoglobin = 6

# Let's assume the BP was recorded as systolic/diastolic, so we'll use 60/40 for clarity.
if systolic_bp > diastolic_bp:
    correct_systolic, correct_diastolic = systolic_bp, diastolic_bp
else:
    correct_systolic, correct_diastolic = diastolic_bp, systolic_bp

assess_patient_and_recommend_treatment(heart_rate, correct_systolic, correct_diastolic, hemoglobin)