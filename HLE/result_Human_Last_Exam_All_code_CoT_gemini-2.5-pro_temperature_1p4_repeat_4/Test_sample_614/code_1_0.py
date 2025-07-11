# Based on the riddle's clues, we deduced the order of the five numbers.
# This script will now print the final sequence.

# Clue 1: "Number 1 is the first..." -> Position 1 is 1
# Clue 5: "...the fifth...lastly follows." -> Position 5 is 5
# Clue 3: "...protects the last;" -> Position 4 is 3
# Clue 4: "...4 ... always goes behind [2]" -> The pair is 2, 4
# This leaves positions 2 and 3 for the pair 2, 4.

# The final sequence is [1, 2, 4, 3, 5]
# Let's print each number of the "final equation" or sequence.

number_1 = 1
number_2 = 2
number_3 = 4
number_4 = 3
number_5 = 5

print("The correct order of the numbers is:")
print(number_1, number_2, number_3, number_4, number_5)