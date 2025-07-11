import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# The exponents derived from physical analysis
n = {
    1: 0,
    2: 1,
    3: -1,
    4: 0,
    5: -2,
    6: -1,
}

# Calculate the sum
total_sum = 0
sum_expression_parts = []
for k in range(1, 7):
    term = k * n[k]
    total_sum += term
    sum_expression_parts.append(f"{k}*({n[k]})")

sum_expression = " + ".join(sum_expression_parts)

# Print the full equation and the final result
print("The values of the exponents are:")
print(f"n1 = {n[1]}, n2 = {n[2]}, n3 = {n[3]}, n4 = {n[4]}, n5 = {n[5]}, n6 = {n[6]}")
print("\nThe final equation is:")
# Manually format for better readability since n_k can be negative
print(f"1*({n[1]}) + 2*({n[2]}) + 3*({n[3]}) + 4*({n[4]}) + 5*({n[5]}) + 6*({n[6]}) = {total_sum}")
print(f"\nThe value of the sum is: {total_sum}")

# Restore original stdout
sys.stdout = original_stdout
# Get the content of the buffer
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)
final_answer = total_sum