# Based on the analysis of the riddle, the sequence of numbers is determined.
# This script will print the final ordered sequence.

# "Number 1 is the first" -> The sequence starts with 1.
num1 = 1

# "Number 4 likes two the most and it always goes behind" -> We have a pair: 2, 4
num2 = 2
num4 = 4

# "Number three is the strongest and it protects the last" -> 3 is before the last number.
num3 = 3

# "the fifth, who never had a thing and lastly follows" -> The last number is 0.
num5 = 0

# Assembling the sequence based on the clues:
# Starts with 1. Ends with 0. 3 is before 0. The pair 2, 4 fills the remaining spots.
# Final sequence: 1, 2, 4, 3, 0
print(f"{num1}, {num2}, {num4}, {num3}, {num5}")