import math

# Define the problem parameters based on the analysis.
num_people = 9
# "black and yellow or blue and white" implies 4 distinct colors.
num_colors = 4

# --- Calculate N (Simultaneous Guessing) ---
# The best strategy guarantees that a number of people equal to the
# floor of the number of people divided by the number of colors will be correct.
# This corresponds to the size of the smallest group in a modulo-based strategy.
N = math.floor(num_people / num_colors)

# --- Calculate M (Sequential Guessing) ---
# When one person speaks first, they can communicate information to the others.
# This strategy ensures all other individuals (num_people - 1) guess correctly.
M = num_people - 1

# --- Calculate the difference M - N ---
difference = M - N

print(f"Maximum guaranteed correct guesses (simultaneous), N: {N}")
print(f"Maximum guaranteed correct guesses (sequential), M: {M}")
print("The difference, representing the additional people who guess correctly, is M - N.")
# The final equation with each number explicitly shown:
print(f"{M} - {N} = {difference}")