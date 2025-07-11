import sys
import io

# capture the output of the print statement
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Start of the actual code ---

# Define the variables based on the scenario
original_estimate = 3000.00
final_invoice = 3500.00

# Under the Ontario Consumer Protection Act, the maximum increase over an estimate is 10%
allowed_increase_decimal = 0.10

# Calculate the maximum allowed increase in dollar value
max_increase_amount = original_estimate * allowed_increase_decimal

# Calculate the total maximum legal amount that can be charged
total_payable_amount = original_estimate + max_increase_amount

print("This script calculates the maximum amount Marc is required to pay under Ontario's Consumer Protection Act.")
print(f"Original estimate: ${original_estimate:,.2f}")
print(f"Final invoice from HR: ${final_invoice:,.2f}")
print(f"Maximum allowed increase over estimate: {allowed_increase_decimal:.0%}")
print("-" * 30)
print(f"The maximum legally chargeable amount is calculated by adding the 10% increase to the original estimate.")
print(f"Therefore, Marc is only required to pay ${total_payable_amount:,.2f}.")
print("\nThe final calculation is:")

# Print the final equation with each number explicitly shown
print(f"{int(original_estimate)} + ({int(original_estimate)} * {allowed_increase_decimal}) = {int(total_payable_amount)}")

# --- End of the actual code ---

# restore the original stdout
sys.stdout = old_stdout
# get the content of the captured output
output = captured_output.getvalue()

# now you can print the output
print(output)