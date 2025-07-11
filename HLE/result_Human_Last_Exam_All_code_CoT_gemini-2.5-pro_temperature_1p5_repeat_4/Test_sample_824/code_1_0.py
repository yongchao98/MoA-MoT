def check_lupus_nephritis_markers(patient_c3, patient_c4):
    """
    Checks patient's complement C3 and C4 levels against normal ranges
    to assess for activity of lupus nephritis.
    """
    # Normal laboratory reference ranges (may vary slightly by lab)
    normal_c3_min = 90  # mg/dL
    normal_c3_max = 180 # mg/dL
    normal_c4_min = 10  # mg/dL
    normal_c4_max = 40  # mg/dL

    print(f"Patient's C3 level: {patient_c3} mg/dL (Normal range: {normal_c3_min}-{normal_c3_max} mg/dL)")
    print(f"Patient's C4 level: {patient_c4} mg/dL (Normal range: {normal_c4_min}-{normal_c4_max} mg/dL)")
    print("-" * 30)

    is_c3_low = patient_c3 < normal_c3_min
    is_c4_low = patient_c4 < normal_c4_min

    if is_c3_low and is_c4_low:
        print("Result: Both C3 and C4 levels are significantly low.")
        print("Interpretation: Low complement levels indicate consumption by the immune system.")
        print("This is a strong indicator of active lupus nephritis, which explains the rapid decline in kidney function.")
    elif is_c3_low:
        print("Result: C3 level is low.")
        print("Interpretation: This suggests complement pathway activation, consistent with a flare of lupus nephritis.")
    elif is_c4_low:
        print("Result: C4 level is low.")
        print("Interpretation: This suggests complement pathway activation, consistent with a flare of lupus nephritis.")
    else:
        print("Result: Complement levels are within the normal range.")
        print("Interpretation: Active, complement-consuming lupus nephritis is less likely.")


# Hypothetical lab values for the patient during her rapid deterioration
patient_c3_value = 40
patient_c4_value = 5

check_lupus_nephritis_markers(patient_c3_value, patient_c4_value)