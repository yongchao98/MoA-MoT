def find_lab_indicator():
    """
    Analyzes a clinical case of rapid renal decline in an SLE patient
    to identify the most indicative lab parameter for its cause.
    """
    # Lab values representing a severe lupus nephritis flare
    patient_c3_level = 35  # Normal is typically 90-180 mg/dL
    patient_c4_level = 8   # Normal is typically 10-40 mg/dL
    patient_creatinine = 9.2 # Normal is < 1.2 mg/dL

    print("Analyzing the clinical scenario to find the lab parameter indicating the CAUSE of rapid renal decline.\n")
    print("The patient's symptoms are characteristic of a severe flare of Systemic Lupus Erythematosus (SLE), leading to rapidly progressive lupus nephritis.")
    print("The direct cause of kidney damage in this condition is massive deposition of immune complexes, which activates the complement system.\n")

    print(f"Patient's Serum Creatinine: {patient_creatinine} mg/dL")
    print("This indicates severe kidney failure. However, it measures the *result* of the damage, not the underlying immunological *cause*.\n")
    
    print("The most direct measure of the ongoing immunological assault is the level of complement proteins, which are consumed in the process.")

    # The "final equation" showing the key parameters
    print("Final Analysis Equation:")
    print("-------------------------")
    print(f"Patient's C3 Level: {patient_c3_level} mg/dL (Severely Decreased)")
    print(f"Patient's C4 Level: {patient_c4_level} mg/dL (Severely Decreased)")
    print("-------------------------")

    print("\nConclusion:")
    print("The combination of severely decreased C3 and C4 levels indicates massive complement consumption.")
    print("This consumption is the hallmark of the immune-complex-mediated damage that caused the patient's rapid progression to end-stage kidney disease.")
    print("\nTherefore, the best lab parameter to indicate the cause of the renal decline is decreased C3 and C4 complement levels.")

find_lab_indicator()
<<<Decreased C3 and C4 complement levels>>>