# Initialize a list to hold the five numbers in order.
# We'll use 0 as a placeholder for unknown numbers.
sequence = [0] * 5

# "Number 1 is the first..."
sequence[0] = 1

# "...the fifth, who... lastly follows."
sequence[4] = 5

# "Number three... protects the last;" (placing it in the 4th position)
sequence[3] = 3

# "Number 4... always goes behind." (implying it follows Number 2)
# The remaining spots are the 2nd and 3rd. To satisfy the clue,
# Number 2 must be in the 2nd position and Number 4 in the 3rd.
sequence[1] = 2
sequence[2] = 4

# Print the final sequence
# The 'sep' argument adds a separator between the printed items.
print(*sequence, sep=', ')
<<<1, 2, 4, 3, 5>>>