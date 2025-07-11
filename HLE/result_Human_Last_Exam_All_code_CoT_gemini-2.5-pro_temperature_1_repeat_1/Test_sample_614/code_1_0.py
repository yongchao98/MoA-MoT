# The riddle describes the order of five specific numbers.
# Based on the clues, the numbers and their positions are determined as follows:

# "Number 1 is the first..." -> The first number is 1.
num1 = 1

# "...the fifth, who never had a thing and lastly follows." -> The last number is 0.
num5 = 0

# "Number three ... protects the last;" -> The fourth number is 3.
num4 = 3

# "Number 4 ... always goes behind [Number 2]."
# This leaves the second and third positions, which must be filled by 2 and 4 respectively.
num2 = 2
num3 = 4

# The final sequence is constructed from these numbers in their determined order.
print(f"{num1} {num2} {num3} {num4} {num5}")