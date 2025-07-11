# The riddle describes the order of the first five natural numbers.
# Based on the clues, the sequence is determined as follows:
# 1. "Number 1 is the first..." -> The first number is 1.
# 2. "...the fifth, who... lastly follows." -> The last number is 5.
# 3. "Number three... protects the last" -> The fourth number is 3.
# 4. "Number 4... always goes behind [two]" -> The numbers 2 and 4 appear as a pair (2, 4).
# This logic leads to the sequence 1, 2, 4, 3, 5.

# The final ordered sequence of numbers.
final_sequence = [1, 2, 4, 3, 5]

# We will print each number in the final sequence, fulfilling the request
# to "output each number in the final equation!".
print("The final order is:")
print(*final_sequence, sep=", ")