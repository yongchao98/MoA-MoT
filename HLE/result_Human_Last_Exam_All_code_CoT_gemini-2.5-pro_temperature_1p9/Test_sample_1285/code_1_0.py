# Part (a) asks for a general formula for the maximum number of roots.
# Based on the analysis of the Wronskian determinant of the given system of functions,
# the maximum number of roots in the interval ]0, 1[ is given by the formula t*(t-1)/2.
part_a_expression = "t*(t-1)/2"

# Part (b) asks for this maximum number when t = 5.
t = 5

# We use the formula from part (a) to calculate the value.
# The final equation is: 5 * (5 - 1) / 2 = result
numerator_val = t * (t - 1)
denominator_val = 2
result_b = numerator_val // denominator_val

# The problem asks to output each number in the final equation.
# We will format the output for part (b) to show the full calculation.
part_b_expression = f"{t}*({t}-1)/{denominator_val} = {result_b}"

# Print the final answer in the specified format: (a) [expression]; (b) [expression].
print(f"(a) {part_a_expression}; (b) {result_b}")