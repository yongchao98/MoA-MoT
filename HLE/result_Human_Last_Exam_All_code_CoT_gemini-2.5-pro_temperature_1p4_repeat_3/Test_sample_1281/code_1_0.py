# Define the initial values from the scenario
estimate_price = 3000.00
invoiced_price = 3500.00
overage_limit_percentage = 10

# Convert the percentage to a decimal for calculation
overage_limit_decimal = overage_limit_percentage / 100.0

# Calculate the maximum allowed price increase in dollars
max_increase_amount = estimate_price * overage_limit_decimal

# Calculate the total maximum legal price
max_legal_price = estimate_price + max_increase_amount

# Print the breakdown of the calculation as per the Ontario Consumer Protection Act
print("This calculation determines the maximum amount Marc is required to pay under Ontario's Consumer Protection Act.")
print("-" * 80)
print(f"Original Estimated Price: ${estimate_price:.2f}")
print(f"Final Invoiced Price: ${invoiced_price:.2f}")
print(f"Under the Act, the price cannot exceed the estimate by more than {overage_limit_percentage}%.")
print("\nStep 1: Calculate the maximum allowed increase.")
print(f"Equation: ${estimate_price:.2f} (Estimate) * {overage_limit_decimal} (Allowed Overage) = ${max_increase_amount:.2f}")

print("\nStep 2: Calculate the total legally enforceable amount.")
print(f"Final Equation: ${estimate_price:.2f} (Estimate) + ${max_increase_amount:.2f} (Max Increase) = ${max_legal_price:.2f}")
print("-" * 80)
print(f"Conclusion: Although the invoice was for ${invoiced_price:.2f}, Marc is only legally required to pay HR ${max_legal_price:.2f}.")
