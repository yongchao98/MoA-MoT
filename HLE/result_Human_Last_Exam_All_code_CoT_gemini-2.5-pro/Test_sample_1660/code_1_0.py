# Define the compound information
tmt_cl_name = "Trimethyltin chloride (TMT-Cl)"
tmt_cl_ld50_oral_rat = 12.6  # mg/kg

tbt_cl_name = "Tributyltin chloride (TBT-Cl)"
tbt_cl_ld50_oral_rat = 129  # mg/kg (using the lower end of the reported range 129-234 mg/kg)

# Explain the concept of LD50
print("The danger or acute toxicity of a chemical is often measured by its LD50 value.")
print("LD50 (Lethal Dose, 50%) is the dose required to be lethal to 50% of a test population.")
print("A lower LD50 value means the substance is more toxic.\n")

# Print the LD50 values for both compounds
print(f"The oral LD50 in rats for {tmt_cl_name} is: {tmt_cl_ld50_oral_rat} mg/kg")
print(f"The oral LD50 in rats for {tbt_cl_name} is: {tbt_cl_ld50_oral_rat} mg/kg\n")

# Compare the values and print the conclusion
toxicity_ratio = tbt_cl_ld50_oral_rat / tmt_cl_ld50_oral_rat

print("To compare their toxicities, we can see how many times larger TBT-Cl's LD50 is than TMT-Cl's LD50.")
print("Calculation:")
print(f"{tbt_cl_ld50_oral_rat} mg/kg (LD50 for TBT-Cl) / {tmt_cl_ld50_oral_rat} mg/kg (LD50 for TMT-Cl) = {toxicity_ratio:.1f}\n")

print(f"This calculation shows that TMT-Cl is approximately {toxicity_ratio:.1f} times more acutely toxic than TBT-Cl.")
print("Therefore, the most important factor is that TMT-Cl has a significantly lower LD50 value, making it much more dangerous.")
