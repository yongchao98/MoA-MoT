# Based on the riddle's clues, the numbers are placed in the following order:
# 1. "Number 1 is the first..." -> 1 is at the start.
# 2. "...the fifth, who ... lastly follows." -> 5 is at the end.
# 3. "Number three ... protects the last;" -> 3 is just before 5.
# 4. "Number 4 ... always goes behind [two]." -> The pair 2, 4 fills the remaining spots.

final_order = [1, 2, 4, 3, 5]

# Print the final sequence of numbers
print("The correct order of the numbers is:")
print(*final_order, sep=", ")