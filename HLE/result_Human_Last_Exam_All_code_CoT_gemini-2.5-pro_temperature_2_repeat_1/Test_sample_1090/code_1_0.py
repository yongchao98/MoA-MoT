def calculate_ldl():
    """
    This function calculates the Low-Density Lipoprotein (LDL) cholesterol
    using the Friedewald equation based on the patient's lipid panel.
    The formula is: LDL = Total Cholesterol - HDL - (Triglycerides / 5).
    This formula is valid for triglyceride levels < 400 mg/dL.
    """
    # Values from the patient's fasting lipid panel
    total_cholesterol = 160  # mg/dL
    hdl = 37  # mg/dL
    triglycerides = 140  # mg/dL

    # Calculate the LDL cholesterol
    # Note: Triglycerides must be divided by 5
    vldl = triglycerides / 5
    calculated_ldl = total_cholesterol - hdl - vldl

    # Print the equation and the result
    print("Calculating LDL using the Friedewald Equation:")
    print(f"LDL = Total Cholesterol - HDL - (Triglycerides / 5)")
    print(f"LDL = {total_cholesterol} - {hdl} - ({triglycerides} / 5)")
    print(f"LDL = {total_cholesterol} - {hdl} - {vldl}")
    print(f"Calculated LDL = {calculated_ldl:.0f} mg/dL")

    # Compare with the reported LDL value from the lab results
    reported_ldl = 95
    print(f"The reported LDL value was {reported_ldl} mg/dL.")
    if int(calculated_ldl) == reported_ldl:
        print("The calculated value matches the reported lab value.")
    else:
        print("The calculated value does not match the reported lab value.")

# Run the calculation
calculate_ldl()
