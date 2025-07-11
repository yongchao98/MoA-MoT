# The Friedewald equation is used to estimate LDL cholesterol:
# LDL = Total Cholesterol - HDL - (Triglycerides / 5)
# All values are in mg/dL.

# Patient's lab values
total_cholesterol = 160
hdl = 37
triglycerides = 140

# Calculate LDL
ldl_c = total_cholesterol - hdl - (triglycerides / 5)

# Print the equation and the result
print("Friedewald Equation for LDL Calculation:")
print(f"{total_cholesterol} (Total Cholesterol) - {hdl} (HDL) - ({triglycerides} (Triglycerides) / 5) = {ldl_c:.0f} mg/dL")

# The calculated LDL of 95 mg/dL matches the patient's lab report.