# Define the financial details from the agreement
purchase_price = 3000
payment_amount = 500
payments_made = 3

# Calculate the total amount Jack has paid
total_paid = payments_made * payment_amount

# In Ontario, the Consumer Protection Act prevents repossession if two-thirds of the price has been paid.
# Let's calculate this threshold.
two_thirds_threshold = (2/3) * purchase_price

print(f"Total Purchase Price: ${purchase_price}")
print(f"Number of payments made: {payments_made}")
print(f"Total Amount Paid by Jack: {payments_made} * ${payment_amount} = ${total_paid}")
print(f"The legal repossession protection threshold (2/3 of the price) is: ${two_thirds_threshold:.2f}")
print("\nComparing the amount paid to the threshold:")

# The final code outputs each number in the final equation as requested.
print(f"Is {total_paid} >= {two_thirds_threshold:.2f}?")
print("The comparison shows that the amount Jack paid is less than the two-thirds threshold.")
print("The final equation is:")
print(total_paid, "<", int(two_thirds_threshold))
print("\nTherefore, Jack does not have statutory protection against repossession based on the amount paid.")
print("However, Gary is not entitled to repossession because he failed to provide proper 'written notice of default' as required by the contract.")
