# Define the LD50 values (oral, mouse) in mg/kg
ld50_tmt_cl = 12.6  # Trimethyltin chloride
ld50_tbt_cl = 129    # Tributyltin chloride

# Print the values for comparison
print("Toxicity Comparison:")
print(f"Trimethyltin chloride (TMT-Cl) LD50 (oral, mouse): {ld50_tmt_cl} mg/kg")
print(f"Tributyltin chloride (TBT-Cl) LD50 (oral, mouse): {ld50_tbt_cl} mg/kg")
print("-" * 30)

# Calculate the toxicity ratio
# A lower LD50 value means higher toxicity.
toxicity_ratio = ld50_tbt_cl / ld50_tmt_cl

# Print the conclusion
print("Conclusion:")
print("A lower LD50 value indicates a higher level of acute toxicity.")
print(f"Based on these values, TMT-Cl is approximately {toxicity_ratio:.1f} times more acutely toxic than TBT-Cl.")
print("This significant difference in LD50 values is the most important quantitative factor showing that TMT-Cl is more dangerous to mammals than TBT-Cl.")
