# The riddle describes the order of five numbers.
# Let's assign the numbers to variables based on the logical deduction.

# "Number 1 is the first..."
num1 = 1

# "...the fifth, who never had a thing and lastly follows."
# The fifth number in the sequence is 0.

# "Number three is the strongest and it protects the last"
# This places 3 just before the last number (0).

# "Number 4 likes two the most and it always goes behind."
# This creates the pair 2, 4.

# Combining the clues, the only possible sequence is 1, 2, 4, 3, 0.
n1 = 1
n2 = 2
n3 = 4
n4 = 3
n5 = 0

# The riddle asks to output each number in the final sequence.
print(f"{n1} {n2} {n3} {n4} {n5}")