# The riddle describes the positions of five numbers.
# Based on the logical deduction from the clues, the final sequence is determined.
# This script stores the solved sequence in a list and prints each number.

solved_sequence = [1, 2, 4, 3, 5]

print("The numbers in the correct order are:")

# We use a loop to print each number from the solved sequence.
# The 'end' parameter is used to print spaces instead of newlines between numbers.
for i, number in enumerate(solved_sequence):
  print(number, end="")
  if i < len(solved_sequence) - 1:
    print(" ", end="")

# Print a final newline for clean output formatting.
print()