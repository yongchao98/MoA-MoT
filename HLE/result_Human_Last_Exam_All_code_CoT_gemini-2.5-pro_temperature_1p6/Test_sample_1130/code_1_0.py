# Player's assets
sestertii = 13
num_vesta_cards = 1

total_houses = 12
brick_houses = 3
num_jupiter_cards = 2

provinces_with_houses = 7
num_saturn_cards = 2

# Goods produced: Brick, Food, Tool, Cloth
num_goods_types_produced = 4
num_mercurius_cards = 2

tool_producing_houses = 3
num_smith_cards = 1 # Smith is a Minerva card for tools

num_colonists = 5
has_concordia_card = True

# --- Score Calculation ---

# 1. Vesta Score (Sestertii)
vesta_points_per_card = sestertii // 10
vesta_score = vesta_points_per_card * num_vesta_cards
print(f"Vesta Score: {vesta_points_per_card} points from sestertii * {num_vesta_cards} Vesta card = {vesta_score} points")

# 2. Jupiter Score (Non-brick cities)
non_brick_cities = total_houses - brick_houses
jupiter_score = non_brick_cities * num_jupiter_cards
print(f"Jupiter Score: {non_brick_cities} non-brick cities * {num_jupiter_cards} Jupiter cards = {jupiter_score} points")

# 3. Saturn Score (Provinces)
saturn_score = provinces_with_houses * num_saturn_cards
print(f"Saturn Score: {provinces_with_houses} provinces * {num_saturn_cards} Saturn cards = {saturn_score} points")

# 4. Mercurius Score (Types of goods)
mercurius_points_per_card = num_goods_types_produced * 2
mercurius_score = mercurius_points_per_card * num_mercurius_cards
print(f"Mercurius Score: ({num_goods_types_produced} types of goods * 2 points) * {num_mercurius_cards} Mercurius cards = {mercurius_score} points")

# 5. Minerva Score (Smith card for Tools)
minerva_points_per_card = tool_producing_houses * 3
minerva_score = minerva_points_per_card * num_smith_cards
print(f"Minerva (Smith) Score: {tool_producing_houses} tool cities * 3 points * {num_smith_cards} Smith card = {minerva_score} points")

# 6. Colonist Score
colonist_score = num_colonists
print(f"Colonist Score: {colonist_score} points from {num_colonists} colonists")

# 7. Concordia Card Score
concordia_score = 7 if has_concordia_card else 0
print(f"Concordia Card Score: {concordia_score} points")

# --- Total Score ---
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + minerva_score + colonist_score + concordia_score

print("\n--- Final Score Calculation ---")
print(f"{vesta_score} (Vesta) + {jupiter_score} (Jupiter) + {saturn_score} (Saturn) + {mercurius_score} (Mercurius) + {minerva_score} (Minerva) + {colonist_score} (Colonists) + {concordia_score} (Concordia) = {total_score}")

print(f"\nThe player's total score is: {total_score}")
<<<70>>>