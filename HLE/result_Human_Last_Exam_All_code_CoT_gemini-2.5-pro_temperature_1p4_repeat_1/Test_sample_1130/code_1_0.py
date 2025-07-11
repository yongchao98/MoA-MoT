# Game state variables
num_houses_brick = 3
num_houses_food = 4
num_houses_tools = 3
num_houses_cloth = 2
num_provinces = 7
sestertii = 13
has_concordia_card = True

# Card counts
num_vesta_cards = 1
num_saturn_cards = 2
num_jupiter_cards = 2
num_mercurius_cards = 2

# --- Score Calculation ---

# 1. Vesta points (1 point per 10 sestertii)
vesta_points_per_card = sestertii // 10
total_vesta_points = num_vesta_cards * vesta_points_per_card

# 2. Saturn points (1 point per province)
saturn_points_per_card = num_provinces
total_saturn_points = num_saturn_cards * saturn_points_per_card

# 3. Jupiter points (1 point per non-brick house)
total_houses = num_houses_brick + num_houses_food + num_houses_tools + num_houses_cloth
non_brick_houses = total_houses - num_houses_brick
jupiter_points_per_card = non_brick_houses
total_jupiter_points = num_jupiter_cards * jupiter_points_per_card

# 4. Mercurius points (2 points per good type)
# The player produces brick, food, tools, and cloth (4 types)
num_good_types = 4
mercurius_points_per_card = num_good_types * 2
total_mercurius_points = num_mercurius_cards * mercurius_points_per_card

# 5. Concordia card points
concordia_points = 7 if has_concordia_card else 0

# 6. Total score
total_score = (total_vesta_points +
               total_saturn_points +
               total_jupiter_points +
               total_mercurius_points +
               concordia_points)

# --- Output ---
print("--- Scoring Breakdown ---")
print(f"Vesta Points: {num_vesta_cards} card * {vesta_points_per_card} points = {total_vesta_points}")
print(f"Saturn Points: {num_saturn_cards} cards * {saturn_points_per_card} points = {total_saturn_points}")
print(f"Jupiter Points: {num_jupiter_cards} cards * {jupiter_points_per_card} points = {total_jupiter_points}")
print(f"Mercurius Points: {num_mercurius_cards} cards * {mercurius_points_per_card} points = {total_mercurius_points}")
print(f"Concordia Card Points: {concordia_points}")
print("---")
print(f"Final Equation: {total_vesta_points} + {total_saturn_points} + {total_jupiter_points} + {total_mercurius_points} + {concordia_points}")
print(f"Total Score: {total_score}")
