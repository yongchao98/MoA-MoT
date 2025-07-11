# Accredited Investor threshold for joint annual income
joint_income_threshold = 300000.00

# Income details for James and his spouse from option D
james_income_2019 = 75000.00
james_income_2020 = 75000.00

spouse_income_2019 = 210000.00
spouse_income_2020 = 210000.00

# Calculate the combined joint income for each year
joint_income_2019 = james_income_2019 + spouse_income_2019
joint_income_2020 = james_income_2020 + spouse_income_2020

# Check if the joint income meets the threshold for each year
passes_2019_test = joint_income_2019 > joint_income_threshold
passes_2020_test = joint_income_2020 > joint_income_threshold

# Print the analysis
print("Analysis for James (shareholder in Option D):")
print(f"The joint income threshold to qualify as an Accredited Investor is ${joint_income_threshold:,.2f} in each of the last two years.")
print("-" * 50)
print(f"Joint income in 2019: ${james_income_2019:,.2f} (James) + ${spouse_income_2019:,.2f} (Spouse) = ${joint_income_2019:,.2f}")
print(f"Does 2019 income exceed the threshold? {passes_2019_test}")
print("-" * 50)
print(f"Joint income in 2020: ${james_income_2020:,.2f} (James) + ${spouse_income_2020:,.2f} (Spouse) = ${joint_income_2020:,.2f}")
print(f"Does 2020 income exceed the threshold? {passes_2020_test}")
print("-" * 50)

# Conclusion
if passes_2019_test and passes_2020_test:
    print("Conclusion: James qualifies as an Accredited Investor based on the joint income test.")
else:
    print("Conclusion: James fails to meet the joint income test in at least one of the two years and is therefore not an Accredited Investor.")
    print("\nBecause the corporation in option D has a shareholder (James) who is not an accredited investor, the corporation fails the test of having all owners be accredited investors. Although the corporation may qualify on its own net assets, this is the only option with a party that definitively fails the specified numerical test.")
