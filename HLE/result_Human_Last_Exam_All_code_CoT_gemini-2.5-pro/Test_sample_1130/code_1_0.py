# Player's assets and cards
sestertii = 13
provinces_with_houses = 7
colonists = 5
total_houses = 12
brick_houses = 3
# The player produces 4 types of goods: brick, food, tool, cloth
num_good_types = 4

# Number of each scoring card
vesta_cards = 1
jupiter_cards = 2
saturn_cards = 2
mercurius_cards = 2
mars_cards = 1  # The "Smith" is a Mars card
concordia_card_points = 7

# --- Score Calculation ---

# 1. Vesta Score (1 VP per 10 sestertii)
vesta_score = vesta_cards * (sestertii // 10)

# 2. Jupiter Score (1 VP per non-brick city)
non_brick_cities = total_houses - brick_houses
jupiter_score = jupiter_cards * non_brick_cities

# 3. Saturnus Score (1 VP per province with a house)
saturnus_score = saturn_cards * provinces_with_houses

# 4. Mercurius Score (2 VP per type of good produced)
mercurius_score = mercurius_cards * num_good_types * 2

# 5. Mars Score (2 VP per colonist)
mars_score = mars_cards * colonists * 2

# 6. Concordia Score (flat 7 VP)
# This is a constant value as the player has the card.

# 7. Total Score
total_score = vesta_score + jupiter_score + saturnus_score + mercurius_score + mars_score + concordia_card_points

# --- Output the results ---
print("Concordia Score Breakdown:")
print(f"Vesta Points: {vesta_score}")
print(f"Jupiter Points: {jupiter_score}")
print(f"Saturnus Points: {saturnus_score}")
print(f"Mercurius Points: {mercurius_score}")
print(f"Mars (Smith) Points: {mars_score}")
print(f"Concordia Card Points: {concordia_card_points}")
print("-" * 20)
# Final equation as requested
print(f"Final Score Calculation: {vesta_score} + {jupiter_score} + {saturnus_score} + {mercurius_score} + {mars_score} + {concordia_card_points} = {total_score}")
print(f"Total Score: {total_score}")
