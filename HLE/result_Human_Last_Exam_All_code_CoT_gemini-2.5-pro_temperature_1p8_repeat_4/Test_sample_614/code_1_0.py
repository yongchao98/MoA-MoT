# Based on the clues from the riddle, we have determined the correct order of the five numbers.
#
# 1. "Number 1 is the first..." -> The sequence starts with 1.
# 2. "...the fifth, who never had a thing and lastly follows." -> The sequence ends with 0.
# 3. "Number three...protects the last" -> 3 comes just before the last number, so ...3, 0.
# 4. "Number 4...always goes behind [two]." -> The pair 2, 4 must appear in that order.
#
# Assembling these pieces gives us the final sequence: 1, 2, 4, 3, 0.
# This script will print each number of that final sequence.

final_sequence = [1, 2, 4, 3, 0]

# Print the final ordered sequence of numbers
print("The final sequence of numbers is:")
print(f"{final_sequence[0]}, {final_sequence[1]}, {final_sequence[2]}, {final_sequence[3]}, {final_sequence[4]}")
