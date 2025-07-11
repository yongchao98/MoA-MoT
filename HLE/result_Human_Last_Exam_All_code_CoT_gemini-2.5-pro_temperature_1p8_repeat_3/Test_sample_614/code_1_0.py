# The riddle describes the order of five natural numbers (1, 2, 3, 4, 5).
# Let's define the final sequence based on the step-by-step analysis of the clues.

# The clues lead to the following order:
# 1st position: Number 1
# 2nd position: Number 2
# 3rd position: Number 4
# 4th position: Number 3
# 5th position: Number 5
final_sequence = [1, 2, 4, 3, 5]

# The sum of the numbers in the sequence.
total = sum(final_sequence)

# The problem asks to output each number in the final equation.
# We will format the sequence as a summation equation.
# We build the string representation of the equation.
# The `map(str, final_sequence)` part converts each number in the list to a string.
# The `' + '.join(...)` part joins them together with ' + ' in between.
equation_str = ' + '.join(map(str, final_sequence))

# Print the final equation.
print(f"{equation_str} = {total}")