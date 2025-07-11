# Player's assets
houses_total = 12
houses_brick = 3
provinces = 7
colonists = 5
sestertii = 13
has_concordia_card = True

# Player's cards for scoring
cards_vesta = 1
cards_jupiter = 2
cards_saturn = 2
cards_mercurius = 2
# The Smith card counts as a Mars card for scoring
cards_mars = 1 

# The types of goods produced are brick, food, tools, and cloth
num_good_types = 4

# --- Score Calculation ---

# 1. Vesta Score (Sestertii)
vesta_score = (sestertii // 10) * cards_vesta
print(f"Vesta Score (Sestertii): ({sestertii} // 10) * {cards_vesta} = {vesta_score}")

# 2. Jupiter Score (Non-brick cities)
non_brick_cities = houses_total - houses_brick
jupiter_score = non_brick_cities * cards_jupiter
print(f"Jupiter Score (Non-brick Cities): ({houses_total} - {houses_brick}) * {cards_jupiter} = {jupiter_score}")

# 3. Saturn Score (Provinces)
saturn_score = provinces * cards_saturn
print(f"Saturn Score (Provinces): {provinces} * {cards_saturn} = {saturn_score}")

# 4. Mercurius Score (Types of Goods)
mercurius_score = num_good_types * cards_mercurius
print(f"Mercurius Score (Goods): {num_good_types} * {cards_mercurius} = {mercurius_score}")

# 5. Mars Score (Colonists)
mars_score = colonists * cards_mars
print(f"Mars Score (Colonists): {colonists} * {cards_mars} = {mars_score}")

# 6. Concordia Card Score
concordia_score = 7 if has_concordia_card else 0
print(f"Concordia Card Score: {concordia_score}")

# --- Total Score ---
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + mars_score + concordia_score
print("\n--- Total Score Calculation ---")
print(f"Total Score: {vesta_score} + {jupiter_score} + {saturn_score} + {mercurius_score} + {mars_score} + {concordia_score} = {total_score}")
<<<53>>>