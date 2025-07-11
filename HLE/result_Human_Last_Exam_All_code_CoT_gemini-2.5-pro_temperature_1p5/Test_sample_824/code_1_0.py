def interpret_lupus_nephritis_labs(patient_c3, patient_c4):
    """
    This function simulates the interpretation of complement levels (C3 and C4)
    in a patient with suspected lupus nephritis flare, explaining their significance.

    Args:
        patient_c3 (int): The patient's measured C3 level in mg/dL.
        patient_c4 (int): The patient's measured C4 level in mg/dL.
    """

    # Typical normal ranges are approximately C3: 90-180 mg/dL and C4: 10-40 mg/dL.
    # We use simplified thresholds for this demonstration.
    c3_low_threshold = 80
    c4_low_threshold = 10

    print("--- Clinical Scenario: Suspected Lupus Nephritis Flare ---")
    print(f"Patient presented with rapidly declining kidney function.")
    print("Investigating for signs of active immune complex deposition.")
    print("\n--- Lab Report: Complement Levels ---")
    print(f"Patient's C3 Level: {patient_c3} mg/dL (Normal range approx. > {c3_low_threshold} mg/dL)")
    print(f"Patient's C4 Level: {patient_c4} mg/dL (Normal range approx. > {c4_low_threshold} mg/dL)")
    print("-" * 35)

    if patient_c3 < c3_low_threshold or patient_c4 < c4_low_threshold:
        print("\n<<< Clinical Interpretation >>>")
        print("Finding: Hypocomplementemia (low C3 and/or C4) is present.")
        print("Conclusion: This finding strongly indicates consumption of complement.")
        print("This consumption is caused by the activation of the classical complement pathway, a process initiated by immune complex deposition in the kidneys.")
        print("Therefore, these decreased levels would have been the best indicator of the underlying cause (active lupus nephritis) for the patient's rapid renal function decline.")
    else:
        print("\n<<< Clinical Interpretation >>>")
        print("Finding: Complement levels are within normal limits.")
        print("Conclusion: Active, complement-consuming immune complex disease is less likely.")
        print("Other causes for renal decline should be prioritized.")


# Simulating the likely lab values for the patient during her acute renal failure.
simulated_patient_c3_level = 45
simulated_patient_c4_level = 6

# Run the interpretation
interpret_lupus_nephritis_labs(simulated_patient_c3_level, simulated_patient_c4_level)