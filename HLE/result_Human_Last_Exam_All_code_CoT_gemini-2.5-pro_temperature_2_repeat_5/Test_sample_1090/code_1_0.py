def calculate_ldl():
    """
    Calculates LDL cholesterol using the Friedewald equation from the patient's lab values.
    The formula is: LDL = Total Cholesterol - HDL - (Triglycerides / 5)
    This calculation is valid when triglycerides are < 400 mg/dL.
    """
    total_cholesterol = 160  # mg/dL
    hdl = 37  # mg/dL
    triglycerides = 140  # mg/dL

    # Calculate the VLDL component (Triglycerides / 5)
    vldl = triglycerides / 5
    
    # Calculate LDL
    ldl = total_cholesterol - hdl - vldl

    # Print the explanation and calculation steps
    print("Calculating LDL Cholesterol using the Friedewald Equation:")
    print(f"LDL = Total Cholesterol - HDL - (Triglycerides / 5)")
    print(f"LDL = {total_cholesterol} - {hdl} - ({triglycerides} / 5)")
    print(f"LDL = {total_cholesterol} - {hdl} - {vldl}")
    print(f"Calculated LDL = {ldl} mg/dL")

# Execute the function
calculate_ldl()