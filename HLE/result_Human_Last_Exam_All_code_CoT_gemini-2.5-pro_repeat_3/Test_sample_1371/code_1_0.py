import math

# The problem constraints lead to a structure of two main blocks, a
# Science/Maths block and an Ethicist/Classicist block. The total number
# of arrangements is calculated by considering the internal permutations
# of these blocks. This involves a case-based analysis based on the genders
# of the people at the interfaces of the blocks.

# The combined coefficient from summing the possibilities of all cases.
# This results from (144 + 48 + 64 + 16) = 272.
coefficient = 272

# The arrangement of the 10 non-rower scientists is a key component. One is at
# the end, and the other 9 can be arranged in 9! ways.
factorial_9 = math.factorial(9)

# The arrangement of the 5 "middle" members of the Ethicist/Classicist block
# (those not at the ends) can be done in 5! ways.
factorial_5 = math.factorial(5)

# The final equation is the product of these numbers.
total_arrangements = coefficient * factorial_9 * factorial_5

# Print the final equation with each number explicitly shown
print("The final equation for the total number of arrangements is:")
print(f"Total Ways = {coefficient} * 9! * 5!")
print(f"Total Ways = {coefficient} * {factorial_9} * {factorial_5}")

# Print the final calculated answer
print("\nThe total number of ways to arrange the table is:")
print(total_arrangements)