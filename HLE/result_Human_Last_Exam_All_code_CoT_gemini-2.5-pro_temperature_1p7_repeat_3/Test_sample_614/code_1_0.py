# This script determines the correct sequence of five numbers based on the provided riddle.
# The core logic is to deduce the identity of the five numbers and then use the clues to place them in the correct order.

# The numbers are identified as {0, 1, 2, 3, 4}.
# The clues about "strength" refer to the length of the English spelling of the numbers.
# We will construct the sequence step by step.

print("Putting the numbers in the correct order based on the riddle's clues:")

# Create a list to hold the final sequence.
final_sequence = ["?"] * 5

# Clue: "the fifth, who never had a thing and lastly follows."
# This sets the last number to 0.
final_sequence[4] = 0

# Clue: "Number three is the strongest and it protects the last;"
# This sets the fourth number to 3.
final_sequence[3] = 3

# Clue: "Number 1 is the first to have something..."
# This sets the first number to 1.
final_sequence[0] = 1

# Clue: "Number 4 likes two the most and it always goes behind."
# This places the remaining numbers, 2 and 4, in the second and third positions respectively.
final_sequence[1] = 2
final_sequence[2] = 4

# The final sequence is complete. Now, we output the result as requested.
# The prompt asks to "output each number in the final equation!".
# We will present the final ordered sequence.

print("The final sequence is:")
print(f"{final_sequence[0]} {final_sequence[1]} {final_sequence[2]} {final_sequence[3]} {final_sequence[4]}")