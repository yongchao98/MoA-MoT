# Define the financial figures from the problem
estimate_price = 3000.00
final_invoice = 3500.00

# The Ontario Consumer Protection Act allows for a maximum 10% increase over the estimate.
allowable_increase_percentage = 0.10

# Calculate the maximum allowable charge under the Ontario CPA
max_allowable_charge = estimate_price * (1 + allowable_increase_percentage)

# Calculate the actual increase amount
actual_increase = final_invoice - estimate_price

print(f"Original Estimated Price: ${estimate_price:.2f}")
print(f"Final Invoiced Price: ${final_invoice:.2f}")
print(f"Actual Price Increase: ${actual_increase:.2f}")
print("\nAccording to the Ontario Consumer Protection Act, the final price cannot exceed the estimate by more than 10%.")
print("Here is the calculation for the maximum allowable charge:")
print(f"${estimate_price:.2f} * (1 + {allowable_increase_percentage:.2f}) = ${max_allowable_charge:.2f}")
