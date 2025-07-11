import math

# --- Player's Game State ---
# Houses and production
houses_total = 12
houses_brick = 3
houses_food = 4
houses_tools = 3
houses_cloth = 2
provinces_with_houses = 7

# Colonists and money
colonists = 5
sestertii = 13

# Storehouse contents
storehouse_cloth = 1
storehouse_tools = 4
storehouse_brick = 1

# Cards
cards_vesta = 1
cards_jupiter = 2
cards_saturn = 2
cards_mercurius = 2
cards_smith = 1 # This is a Minerva card for Tools
has_concordia_card = True

# --- Game Rules Data ---
good_values = {
    "brick": 3,
    "food": 4,
    "tools": 5,
    "cloth": 6,
    "wine": 7
}

total_score = 0

print("Calculating Concordia Score:\n")

# 1. Vesta Score
vesta_goods_value = (storehouse_cloth * good_values["cloth"]) + \
                    (storehouse_tools * good_values["tools"]) + \
                    (storehouse_brick * good_values["brick"])
vesta_total_wealth = vesta_goods_value + sestertii
vesta_points_per_card = math.floor(vesta_total_wealth / 10)
vesta_score = cards_vesta * vesta_points_per_card
total_score += vesta_score
print(f"Vesta Points: {cards_vesta} card * floor(({storehouse_cloth} cloth * {good_values['cloth']} + {storehouse_tools} tools * {good_values['tools']} + {storehouse_brick} brick * {good_values['brick']} + {sestertii} sestertii) / 10) = {vesta_score}")


# 2. Jupiter Score
non_brick_houses = houses_total - houses_brick
jupiter_score = cards_jupiter * non_brick_houses
total_score += jupiter_score
print(f"Jupiter Points: {cards_jupiter} cards * ({houses_total} total houses - {houses_brick} brick houses) = {jupiter_score}")

# 3. Saturn Score
saturn_provinces_sets = math.floor(provinces_with_houses / 2)
saturn_score = cards_saturn * saturn_provinces_sets
total_score += saturn_score
print(f"Saturn Points: {cards_saturn} cards * floor({provinces_with_houses} provinces / 2) = {saturn_score}")

# 4. Mercurius Score
produced_good_types = 0
if houses_brick > 0: produced_good_types += 1
if houses_food > 0: produced_good_types += 1
if houses_tools > 0: produced_good_types += 1
if houses_cloth > 0: produced_good_types += 1
# Assuming no wine houses were mentioned
mercurius_points_per_card = produced_good_types * 2
mercurius_score = cards_mercurius * mercurius_points_per_card
total_score += mercurius_score
print(f"Mercurius Points: {cards_mercurius} cards * {produced_good_types} types of goods * 2 = {mercurius_score}")

# 5. Minerva Score (Smith)
smith_points = cards_smith * houses_tools * good_values["tools"]
total_score += smith_points
print(f"Minerva (Smith) Points: {cards_smith} card * {houses_tools} tool cities * {good_values['tools']} points per city = {smith_points}")

# 6. Concordia Card Score
concordia_score = 7 if has_concordia_card else 0
total_score += concordia_score
print(f"Concordia Card Points: {concordia_score}")


print("\n--------------------")
print(f"Final Score Equation: {vesta_score} + {jupiter_score} + {saturn_score} + {mercurius_score} + {smith_points} + {concordia_score}")
print(f"Total Score: {total_score}")
print("--------------------")

<<<66>>>