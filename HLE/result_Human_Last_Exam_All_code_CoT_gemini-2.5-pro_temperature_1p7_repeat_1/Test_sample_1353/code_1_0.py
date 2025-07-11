# Define the initial bi-degree for part (a)
a_start = 4
b_start = 3

# Calculate the string parameter k
k = a_start - b_start

# Calculate the bi-degree of the terminal polynomial
a_end = a_start - k
b_end = b_start + k

# Format the answers for all parts.
# For part (a), show the full calculation in the final expression as requested.
answer_a_expression = f"({a_start} - {k}, {b_start} + {k}) = ({a_end}, {b_end})"
answer_b_expression = "sum_{i=1 to b} r_i = b"
answer_c_expression = "Yes"

# Print the final result in the specified format
print(f"a) {answer_a_expression} b) {answer_b_expression} c) {answer_c_expression}")