def calculate_ldl():
    """
    Calculates LDL cholesterol using the Friedewald equation based on the patient's lab values.
    The Friedewald equation is: LDL = Total Cholesterol - HDL - (Triglycerides / 5).
    This calculation is generally accurate when triglyceride levels are below 400 mg/dL.
    """
    # Patient's fasting lipid panel values from the case study
    total_cholesterol = 160  # mg/dL
    hdl = 37  # mg/dL
    triglycerides = 140  # mg/dL

    # Perform the calculation
    ldl_calculated = total_cholesterol - hdl - (triglycerides / 5)

    # Output the explanation, the equation with the patient's values, and the result.
    print("The Friedewald equation is used to estimate LDL cholesterol from a lipid panel.")
    print("Formula: LDL = Total Cholesterol - HDL - (Triglycerides / 5)")
    print("Based on the patient's lab results, the calculation is:")
    print(f"LDL = {total_cholesterol} - {hdl} - ({triglycerides} / 5)")
    print(f"Calculated LDL Cholesterol = {ldl_calculated} mg/dL")
    print("\nThis calculated value matches the reported LDL of 95 mg/dL in the patient's chart.")

if __name__ == '__main__':
    calculate_ldl()