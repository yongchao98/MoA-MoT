# Define the financial figures from the scenario
estimate_price = 3000.00
invoiced_price = 3500.00

# Define the relevant legal parameter from Ontario's Consumer Protection Act
cpa_allowable_overage_rate = 0.10  # 10%

# Explain the legal principle
print("According to Ontario's Consumer Protection Act, the final price charged to a consumer cannot exceed the written estimate by more than 10%.")
print("First, let's calculate the maximum allowable increase based on the estimate.")

# Calculate the maximum dollar value of the increase
max_increase_amount = estimate_price * cpa_allowable_overage_rate

# Calculate the total maximum price that can be legally charged
max_legal_price = estimate_price + max_increase_amount

# Display the step-by-step calculation
print(f"\nOriginal Estimate: ${estimate_price:.2f}")
print(f"Maximum Allowable Increase (10%): ${max_increase_amount:.2f}")
print(f"Total Maximum Legal Price: ${max_legal_price:.2f}")

# Display the final equation showing all numbers as required
print("\nThe final equation for the amount Marc is required to pay is:")
print(f"${estimate_price:.2f} (Estimate) + (${estimate_price:.2f} * {cpa_allowable_overage_rate}) = ${max_legal_price:.2f}")

# State the final conclusion
print(f"\nSince the invoice for ${invoiced_price:.2f} is higher than the maximum legal price of ${max_legal_price:.2f}, Marc is only required to pay ${max_legal_price:.2f}.")
print("This reasoning directly supports Answer B.")
