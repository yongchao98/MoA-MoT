import math

# Player's assets and cards based on the scenario
houses_brick = 3
houses_food = 4
houses_tools = 3
houses_cloth = 2
provinces_with_houses = 7
sestertii = 13
has_concordia_card = True

# Cards in player's hand
cards_vesta = 1
cards_jupiter = 2
cards_saturn = 2
cards_mercurius = 2
cards_smith = 1 # The Smith is a Minerva card for tools

print("--- Concordia Final Score Calculation ---")

# 1. Vesta Score (from Sestertii)
# Vesta cards grant 1 point for every 10 sestertii.
vesta_factor = sestertii // 10
vesta_score = cards_vesta * vesta_factor
print(f"\nVesta Score (for money):")
print(f"Calculation: {cards_vesta} Vesta card * ({sestertii} sestertii // 10) = {vesta_score} points")

# 2. Jupiter Score (from non-brick cities)
# Jupiter cards grant 1 point per city that does NOT produce brick.
non_brick_cities = houses_food + houses_tools + houses_cloth
jupiter_score = cards_jupiter * non_brick_cities
print(f"\nJupiter Score (for non-brick cities):")
print(f"Calculation: {cards_jupiter} Jupiter cards * {non_brick_cities} cities ({houses_food} Food + {houses_tools} Tools + {houses_cloth} Cloth) = {jupiter_score} points")

# 3. Saturn Score (from provinces)
# Saturn cards grant 1 point per province with at least one house.
saturn_score = cards_saturn * provinces_with_houses
print(f"\nSaturn Score (for provinces):")
print(f"Calculation: {cards_saturn} Saturn cards * {provinces_with_houses} provinces = {saturn_score} points")

# 4. Mercurius Score (for types of goods)
# Mercurius cards grant 2 points for each type of good produced.
# The player produces brick, food, tools, and cloth, which is 4 types.
types_of_goods = 4
mercurius_score = cards_mercurius * types_of_goods * 2
print(f"\nMercurius Score (for types of goods):")
print(f"Calculation: {cards_mercurius} Mercurius cards * {types_of_goods} types of goods * 2 points/type = {mercurius_score} points")

# 5. Minerva (Smith) Score (from specialist cities)
# The Smith card grants 3 points per tool-producing city.
smith_score = cards_smith * houses_tools * 3
print(f"\nSmith Score (for Tool cities):")
print(f"Calculation: {cards_smith} Smith card * {houses_tools} Tool cities * 3 points/city = {smith_score} points")

# 6. Concordia Card Score
# The Concordia card gives a flat 7 points.
concordia_score = 7 if has_concordia_card else 0
print(f"\nConcordia Card Score:")
print(f"Calculation: The card is worth a flat {concordia_score} points.")

# 7. Final Total
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + smith_score + concordia_score

print("\n-------------------------------------------")
print("           Final Score Summary")
print("-------------------------------------------")
print(f"Total Score = Vesta + Jupiter + Saturn + Mercurius + Smith + Concordia")
print(f"Total Score = {vesta_score} + {jupiter_score} + {saturn_score} + {mercurius_score} + {smith_score} + {concordia_score}")
print(f"Final total = {total_score} points")

<<<65>>>