import decimal

# Define the financial values from the problem
estimate_price = decimal.Decimal('3000.00')
invoiced_price = decimal.Decimal('3500.00')
cpa_limit_percentage = decimal.Decimal('0.10') # 10% limit

# Explain the legal reasoning
print("Based on Ontario's Consumer Protection Act (CPA), a business cannot charge more than 10% above its written estimate.")
print("While the work was done in Quebec, the contract was made with an Ontario-based company in Ontario, making the company subject to the CPA's rules on estimates.")
print("-" * 50)

# Calculate the maximum legally chargeable amount
allowed_increase_amount = estimate_price * cpa_limit_percentage
max_legal_price = estimate_price + allowed_increase_amount

# Display the step-by-step calculation
print("Here is the calculation based on the 10% rule:")
print(f"Original Estimate: ${estimate_price:.2f}")
print(f"Maximum Allowed Increase (10% of estimate): ${allowed_increase_amount:.2f}")

# Display the final equation with all its numbers
# We use .to_integral_value() for cleaner presentation in the final equation string
est_int = estimate_price.to_integral_value()
inc_int = allowed_increase_amount.to_integral_value()
max_int = max_legal_price.to_integral_value()

print(f"\nFinal Equation: ${est_int} + ${inc_int} = ${max_int}")

# State the final advice for Marc
print("-" * 50)
print(f"The invoiced price of ${invoiced_price:.2f} is higher than the legally permissible amount.")
print(f"Therefore, Marc is only required to pay HR ${max_legal_price:.2f}.")
