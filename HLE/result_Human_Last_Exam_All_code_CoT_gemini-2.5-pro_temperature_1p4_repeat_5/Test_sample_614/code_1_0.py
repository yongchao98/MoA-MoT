# Plan:
# 1. Define the numbers in the final sequence based on the riddle's clues.
#    - Number 1 is first.
#    - Number 5 is last.
#    - Number 3 is before the last.
#    - Number 2 is followed by Number 4.
#    - This gives the sequence: 1, 2, 4, 3, 5
# 2. The prompt asks to "output each number in the final equation!".
#    I will interpret this as printing the final sequence in a clear format.

pos1 = 1
pos2 = 2
pos3 = 4
pos4 = 3
pos5 = 5

# Print the final sequence, showing each number.
print("The final sequence of the numbers is an equation of order:")
print(f"{pos1} -> {pos2} -> {pos3} -> {pos4} -> {pos5}")
