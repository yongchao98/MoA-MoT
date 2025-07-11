# The transaction is governed by the Ontario Consumer Protection Act because the business is located in Ontario.
# The Act states that an invoice cannot exceed a written estimate by more than 10%.

# Define the initial values from the scenario
original_estimate = 3000.00
invoiced_amount = 3500.00
overage_limit_percentage = 0.10 # This represents the 10% limit

# Calculate the maximum allowable amount over the estimate
max_allowable_overage = original_estimate * overage_limit_percentage

# Calculate the maximum legal price that can be invoiced
max_legal_price = original_estimate + max_allowable_overage

# Print the step-by-step breakdown of the calculation
print(f"Original Estimate: ${original_estimate:.2f}")
print(f"Invoiced Amount: ${invoiced_amount:.2f}")
print("---")
print("Under the Ontario Consumer Protection Act, an invoice cannot exceed an estimate by more than 10%.")
print("Calculating the maximum legal charge:")
print(f"Maximum allowed increase = {int(overage_limit_percentage * 100)}% of the estimate")
print(f"Equation for the legally required payment: ${original_estimate:.2f} (Estimate) + (${original_estimate:.2f} * {overage_limit_percentage}) = ${max_legal_price:.2f}")
print("---")
print(f"Since the invoiced amount of ${invoiced_amount:.2f} is greater than the maximum legal charge of ${max_legal_price:.2f}, Marc is only required to pay HR ${max_legal_price:.2f}.")
