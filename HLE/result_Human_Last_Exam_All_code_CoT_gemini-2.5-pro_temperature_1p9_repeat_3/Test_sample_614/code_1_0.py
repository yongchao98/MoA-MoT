# The riddle is a logic puzzle. Here is the step-by-step deduction encoded.

# Start with 5 unknown positions.
sequence = [None, None, None, None, None]

# Clue: "...the fifth, who ... lastly follows."
# This places 5 at the last position (index 4).
sequence[4] = 5

# Clue: "Number three ... protects the last;"
# This places 3 just before the last element (index 3).
sequence[3] = 3

# Clue: "Number 2 is the best of them;" interpreted as #2 being first.
# This places 2 at the first position (index 0).
sequence[0] = 2

# Clue: "Number 4 ... always goes behind [2]."
# This places 4 immediately after 2 (index 1).
sequence[1] = 4

# The only remaining number is 1, and the only open spot is at index 2.
sequence[2] = 1

# The final resolved sequence.
# Let's print each number in the final equation.
print(*sequence)