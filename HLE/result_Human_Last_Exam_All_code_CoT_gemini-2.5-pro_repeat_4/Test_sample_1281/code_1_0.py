# Define the known financial figures from the scenario.
estimate_price = 3000.00
invoiced_price = 3500.00

# The Ontario Consumer Protection Act limits price increases over an estimate to 10%.
cpa_limit_percentage = 0.10

# Calculate the maximum allowable increase in dollars.
allowable_increase = estimate_price * cpa_limit_percentage

# Calculate the maximum legal price that can be charged.
max_legal_price = estimate_price + allowable_increase

# Print the explanation and the final calculation.
print("Under Ontario's Consumer Protection Act, a final price cannot exceed a given estimate by more than 10%.")
print(f"Original Estimate: ${estimate_price:.2f}")
print(f"Maximum Allowable Increase (10%): ${allowable_increase:.2f}")
print("\nCalculating the maximum legal price:")
print(f"${estimate_price:.2f} (Estimate) + ${allowable_increase:.2f} (10% Increase) = ${max_legal_price:.2f} (Maximum Legal Price)")
print(f"\nSince the invoice was for ${invoiced_price:.2f}, Marc is only legally required to pay ${max_legal_price:.2f}.")
