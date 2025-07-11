# Define the initial values from the problem
original_estimate = 3000.00
invoiced_amount = 3500.00
allowable_overage_percentage = 0.10 # 10% as per the Consumer Protection Act

# Calculate the maximum allowable increase in dollars
max_increase_amount = original_estimate * allowable_overage_percentage

# Calculate the total maximum amount Marc is required to pay
max_payable_amount = original_estimate + max_increase_amount

# Explain the rule and show the calculation
print("Under Ontario's Consumer Protection Act, a price cannot exceed an estimate by more than 10%.")
print("The calculation for the maximum payable amount is as follows:")
print(f"Original Estimate: ${original_estimate:,.2f}")
print(f"Allowable Increase (10%): ${max_increase_amount:,.2f}")
print("\nFinal Equation:")
print(f"${int(original_estimate)} + ${int(max_increase_amount)} = ${int(max_payable_amount)}")

print(f"\nTherefore, Marc is required to pay HR ${max_payable_amount:,.2f}, not the invoiced amount of ${invoiced_amount:,.2f}.")