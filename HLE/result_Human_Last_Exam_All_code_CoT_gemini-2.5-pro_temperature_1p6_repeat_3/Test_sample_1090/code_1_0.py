def calculate_ldl_from_panel():
    """
    Calculates LDL cholesterol using the Friedewald equation from the patient's lab values.
    The equation is valid for Triglycerides < 400 mg/dL.
    """
    # Patient's lipid panel values from the case description
    total_cholesterol = 160  # mg/dL
    hdl = 37  # mg/dL
    triglycerides = 140  # mg/dL

    # The Friedewald equation estimates VLDL as Triglycerides / 5
    vldl_estimate = triglycerides / 5

    # Calculate LDL
    calculated_ldl = total_cholesterol - hdl - vldl_estimate

    # Print the equation with the patient's values
    print("Calculating LDL Cholesterol using the Friedewald Equation:")
    print(f"LDL = Total Cholesterol ({total_cholesterol}) - HDL ({hdl}) - (Triglycerides ({triglycerides}) / 5)")
    
    # Print the final calculation and result
    print(f"LDL = {total_cholesterol} - {hdl} - {int(vldl_estimate)}")
    print(f"Calculated LDL = {int(calculated_ldl)} mg/dL")
    print("(This matches the reported LDL of 95 mg/dL)")

calculate_ldl_from_panel()