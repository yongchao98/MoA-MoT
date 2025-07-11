def assess_trauma_patient():
    """
    Analyzes patient data to diagnose shock and recommend first-line treatment.
    """
    # Patient's vital signs and lab values
    heart_rate = 160  # beats per minute
    # Blood pressure is likely 60/40 mmHg; using 60 for systolic.
    systolic_bp = 60  # mmHg
    hemoglobin = 6  # gm/dL
    clinical_signs = ["profuse bleeding", "cold and clammy skin", "disorientation"]

    print("Patient's Key Data:")
    print(f"Heart Rate: {heart_rate} bpm")
    print(f"Systolic Blood Pressure: {systolic_bp} mmHg")
    print(f"Hemoglobin: {hemoglobin} gm/dL")
    print(f"Clinical Signs: {', '.join(clinical_signs)}\n")

    # Calculate Shock Index (SI) = Heart Rate / Systolic BP
    # A normal SI is 0.5-0.7. An SI > 1.0 indicates severe shock.
    shock_index = heart_rate / systolic_bp
    
    print("Calculating Shock Index (Heart Rate / Systolic BP):")
    print(f"Equation: {heart_rate} / {systolic_bp} = {shock_index:.2f}\n")

    # Diagnosis based on data
    if shock_index > 1.0 and systolic_bp < 90 and hemoglobin < 7:
        print("Diagnosis: The patient is in severe hemorrhagic shock due to massive blood loss.")
        print("The extremely high Shock Index confirms a state of circulatory collapse requiring immediate intervention.\n")
    
    # Recommended Treatment
    print("First-Line Treatment Recommendation:")
    print("The immediate priority is to restore circulating blood volume to improve tissue perfusion.")
    print("The standard of care is rapid intravenous resuscitation with isotonic crystalloids.\n")
    print("Answer Choice C: Intravenous resuscitation of normal saline or Ringer's lactate")

assess_trauma_patient()