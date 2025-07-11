def calculate_ldl():
    """
    Calculates LDL cholesterol using the Friedewald equation.
    The formula is: LDL = Total Cholesterol - HDL - (Triglycerides / 5)
    This formula is generally valid if the triglyceride level is less than 400 mg/dL.
    """
    total_cholesterol = 160  # mg/dL
    hdl = 37  # mg/dL
    triglycerides = 140  # mg/dL

    # Check if the Friedewald equation is applicable
    if triglycerides >= 400:
        print("Triglyceride level is too high for an accurate Friedewald calculation.")
        return

    # Calculate LDL
    ldl_calculated = total_cholesterol - hdl - (triglycerides / 5)

    # Print the equation and the result
    print("Verifying the LDL cholesterol using the Friedewald equation:")
    print(f"LDL = Total Cholesterol - HDL - (Triglycerides / 5)")
    print(f"LDL = {total_cholesterol} - {hdl} - ({triglycerides} / 5)")
    print(f"Calculated LDL = {int(ldl_calculated)} mg/dL")
    print(f"This matches the reported LDL of 95 mg/dL, indicating his lipid-lowering therapy is effective.")

calculate_ldl()