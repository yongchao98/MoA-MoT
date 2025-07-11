# This script determines the correct order of numbers based on the riddle's clues.
# The logic is as follows:
# 1. The number 5 is last ("lastly follows").
# 2. The number 3 is just before the last number ("protects the last").
# 3. The numbers 2 and 4 appear as a pair, 2 followed by 4 ("4...goes behind [2]").
#    This pair must fill the first two available slots.
# 4. The number 1 is the only number left and fills the only remaining slot.

# This results in the sequence [2, 4, 1, 3, 5].
# The script will now print this final sequence.

final_sequence = [2, 4, 1, 3, 5]

# We are asked to output each number.
# Let's print them separated by spaces.
print("The final ordered sequence is:")
print(*final_sequence)