# This script calculates the maximal possible score difference (X-Y) in a single Tichu round
# under the constraint that the winning team does not finish 1st and 2nd.

# --- Step 1: Define the components for the winning team's score (X) ---

# To maximize the difference, the winning team (Team A) must win all 100 card points.
x_cards = 100

# To maximize the bonus, one player on Team A makes a successful "Grand Tichu" call (+200),
# while their partner makes no call. This requires the calling player to finish first.
x_bonus = 200

# Calculate the total score for the winning team (X).
X = x_cards + x_bonus

# --- Step 2: Define the components for the losing team's score (Y) ---

# As Team A won all card points, the losing team (Team B) has 0 card points.
y_cards = 0

# To minimize the losing team's score, both players on Team B make "Grand Tichu" calls and fail.
# This results in a -200 penalty for each player.
y_bonus = -200 + -200

# Calculate the total score for the losing team (Y).
Y = y_cards + y_bonus

# --- Step 3: Calculate the maximal difference X - Y ---

difference = X - Y

# --- Step 4: Print the detailed calculation and the final answer ---

print("Scenario for Maximal Score Difference:")
print(f"Winning Team's Score (X) = Card Points + Bonus Points")
print(f"X = {x_cards} + {x_bonus} = {X}")
print("")
print(f"Losing Team's Score (Y) = Card Points + Bonus Points")
print(f"Y = {y_cards} + ({y_bonus}) = {Y}")
print("")
print("Final Calculation of the Maximal Difference (X - Y):")
print(f"{X} - ({Y}) = {difference}")
