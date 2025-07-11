def calculate_ldl():
    """
    Calculates LDL cholesterol using the Friedewald equation based on the patient's lab values.
    The Friedewald equation is: LDL = Total Cholesterol - HDL - (Triglycerides / 5).
    This formula is generally valid for Triglyceride levels < 400 mg/dL.
    """
    total_cholesterol = 160  # mg/dL
    hdl = 37  # mg/dL
    triglycerides = 140  # mg/dL

    # Check if the calculation is valid
    if triglycerides >= 400:
        print("Triglyceride level is too high for an accurate Friedewald calculation.")
        return

    # Calculate VLDL (Very Low-Density Lipoprotein)
    vldl = triglycerides / 5
    
    # Calculate LDL (Low-Density Lipoprotein)
    ldl_calculated = total_cholesterol - hdl - vldl

    print("Calculating LDL using the Friedewald equation: LDL = Total Cholesterol - HDL - (Triglycerides / 5)")
    print(f"LDL = {total_cholesterol} - {hdl} - ({triglycerides} / 5)")
    print(f"Calculated LDL Cholesterol: {int(ldl_calculated)} mg/dL")

# Execute the function to show the calculation
calculate_ldl()