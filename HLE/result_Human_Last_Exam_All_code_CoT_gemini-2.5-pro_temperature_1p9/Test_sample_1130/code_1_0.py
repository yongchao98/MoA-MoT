# Player's assets
sestertii = 13
total_houses = 12
brick_houses = 3
provinces_with_houses = 7
goods_produced_types = 4  # Brick, Food, Tools, Cloth
colonists = 5

# Player's cards
vesta_cards = 1
jupiter_cards = 2
saturn_cards = 2
mercurius_cards = 2
smith_cards = 1  # Scores like a Mars card
concordia_card = 1

# 1. Vesta Score (1 VP per 10 sestertii)
vesta_points_per_card = sestertii // 10
vesta_score = vesta_cards * vesta_points_per_card
print(f"Vesta Score: {vesta_cards} card * ({sestertii} sestertii // 10) = {vesta_score}")

# 2. Jupiter Score (1 VP per non-brick city)
non_brick_cities = total_houses - brick_houses
jupiter_score = jupiter_cards * non_brick_cities
print(f"Jupiter Score: {jupiter_cards} cards * ({total_houses} total houses - {brick_houses} brick houses) = {jupiter_score}")

# 3. Saturnus Score (1 VP per province)
saturn_score = saturn_cards * provinces_with_houses
print(f"Saturnus Score: {saturn_cards} cards * {provinces_with_houses} provinces = {saturn_score}")

# 4. Mercurius Score (2 VP per type of good)
vp_per_good_type = 2
mercurius_score = mercurius_cards * goods_produced_types * vp_per_good_type
print(f"Mercurius Score: {mercurius_cards} cards * {goods_produced_types} types of goods * {vp_per_good_type} VP = {mercurius_score}")

# 5. Smith (Mars) Score (2 VP per colonist)
vp_per_colonist = 2
smith_score = smith_cards * colonists * vp_per_colonist
print(f"Smith (Mars) Score: {smith_cards} card * {colonists} colonists * {vp_per_colonist} VP = {smith_score}")

# 6. Concordia Card Score
concordia_score = 7 * concordia_card
print(f"Concordia Card Score: {concordia_score}")

# 7. Total Score
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + smith_score + concordia_score
print("\n--- Final Score Calculation ---")
print(f"Total Score = {vesta_score} (Vesta) + {jupiter_score} (Jupiter) + {saturn_score} (Saturnus) + {mercurius_score} (Mercurius) + {smith_score} (Smith) + {concordia_score} (Concordia)")
print(f"Final Score: {total_score}")

<<<66>>>