# The final sequence is determined by the clues in the riddle.
# 1. The number 5 is last.
# 2. The number 3 "protects the last", so it is second-to-last.
# 3. The number 1 is "the first".
# 4. The number 4 "goes behind" the number 2, filling the remaining spots.
# This results in the sequence: 1, 2, 4, 3, 5.

final_sequence = [1, 2, 4, 3, 5]

# As requested, here is the code to output each number in the final sequence.
print("The correct order is:")

# We use a loop to print each number from the list.
# The 'end=" "' argument prints a space instead of a newline after each number.
for number in final_sequence:
    print(number, end=" ")

# Print a final newline character for clean formatting.
print()