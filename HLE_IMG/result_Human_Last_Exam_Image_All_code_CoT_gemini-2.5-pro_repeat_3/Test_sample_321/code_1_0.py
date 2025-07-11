import collections

# Mapping of plot number to the identified parameter identifier
# p_n: identifier of the parameter varied in plot n
# Plot 1: a_i -> 5
# Plot 2: q_f -> 15
# Plot 3: q_s -> 13
# Plot 4: q_l -> 14
# Plot 5: f_s -> 6
# Plot 6: beta_h -> 9
# Plot 7: c_l -> 7
# Plot 8: mu_h -> 8
# Plot 9: mu_n -> 3
solution_map = {
    1: 5,
    2: 15,
    3: 13,
    4: 14,
    5: 6,
    6: 9,
    7: 7,
    8: 8,
    9: 3
}

# The problem asks for the answer as a sequence {p_1, p_2, ..., p_9}
# We will create an ordered list of the identifiers.
ordered_solution = [solution_map[i] for i in range(1, 10)]

# Print the final sequence as requested
# The problem mentions "final equation", which we interpret as the final sequence of numbers.
final_answer_string = "{" + ", ".join(map(str, ordered_solution)) + "}"
print(f"The final sequence of parameter identifiers is:")
print(final_answer_string)