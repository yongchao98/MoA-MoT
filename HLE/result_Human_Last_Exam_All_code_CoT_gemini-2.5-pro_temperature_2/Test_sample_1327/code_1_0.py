# The given sequence
sequence = [2, 11, 23, 51, 119]

# We observe the pattern x_n+1 = 2 * x_n + a_n starting from the second term.
# Let's find the sequence a_n.
a_n = []
# Start from the second term of the sequence (index 1)
for i in range(1, len(sequence) - 1):
  a = sequence[i+1] - 2 * sequence[i]
  a_n.append(a)

# The sequence a_n is [1, 5, 17]

# Now, let's find the pattern in a_n by looking at the differences.
differences = []
for i in range(len(a_n) - 1):
  diff = a_n[i+1] - a_n[i]
  differences.append(diff)
  
# The differences are [4, 12]

# This is a geometric progression with a ratio of 3.
ratio = differences[1] / differences[0]

# Calculate the next difference.
next_difference = differences[-1] * ratio

# Calculate the next term in the a_n sequence.
next_a = a_n[-1] + next_difference

# Get the last term of the original sequence.
last_term = sequence[-1]

# Calculate the final result.
result = 2 * last_term + next_a

# Print the final equation with all its components.
# The formula is 2 * last_term + next_a = result
print(f"2 * {last_term} + {int(next_a)} = {int(result)}")