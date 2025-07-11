def calculate_and_explain_ldl():
    """
    This function calculates the Low-Density Lipoprotein (LDL) Cholesterol
    using the Friedewald equation from the patient's lipid panel.
    The formula is: LDL = Total Cholesterol - HDL - (Triglycerides / 5)
    This is done to fulfill the request for a coded calculation based on the case data.
    """
    total_cholesterol = 160
    hdl = 37
    triglycerides = 140

    # Calculate the LDL value
    calculated_ldl = total_cholesterol - hdl - (triglycerides / 5)

    print("Calculating LDL Cholesterol using the Friedewald Equation:")
    print(f"LDL = Total Cholesterol - HDL - (Triglycerides / 5)")
    print(f"LDL = {total_cholesterol} - {hdl} - ({triglycerides} / 5)")
    print(f"LDL = {total_cholesterol} - {hdl} - {int(triglycerides/5)}")
    print(f"Calculated LDL is: {int(calculated_ldl)} mg/dL")
    print("\nThis calculated value matches the provided lab result of 95 mg/dL.")
    print("\nHowever, the primary clinical issue is the patient's failure to progress in physical therapy, which is most fundamentally impacted by his poor nutritional status (BMI 18.5). Therefore, the most appropriate next step is a nutritional consultation.")

calculate_and_explain_ldl()