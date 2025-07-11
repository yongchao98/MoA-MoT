import math

# Player's assets
houses_brick = 3
houses_food = 4
houses_tools = 3
houses_cloth = 2
total_houses = houses_brick + houses_food + houses_cloth + houses_tools
provinces_with_houses = 7
sestertii = 13
has_concordia_card = True

# Storehouse contents
storehouse_cloth = 1
storehouse_tools = 4
storehouse_brick = 1

# Player's cards
cards_saturn = 2
cards_jupiter = 2
cards_vesta = 1
cards_mercurius = 2
cards_smith = 1 # Minerva card for tools

# --- Score Calculation ---

# 1. Vesta Score
goods_value = (storehouse_cloth * 6) + (storehouse_tools * 5) + (storehouse_brick * 3)
vesta_base_value = goods_value + sestertii
vesta_score = math.floor(vesta_base_value / 10) * cards_vesta

# 2. Jupiter Score
non_brick_houses = total_houses - houses_brick
jupiter_score = non_brick_houses * cards_jupiter

# 3. Saturn Score
saturn_score = provinces_with_houses * cards_saturn

# 4. Mercurius Score
produced_goods_types = 4 # brick, food, tools, cloth
mercurius_score = (produced_goods_types * 2) * cards_mercurius

# 5. Smith (Minerva) Score
# Smith card gives 3 points per house in a tool-producing city
smith_score = houses_tools * 3 * cards_smith

# 6. Concordia Card Score
concordia_score = 7 if has_concordia_card else 0

# 7. Total Score
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + smith_score + concordia_score

# --- Output the final equation ---
print(f"Vesta ({vesta_score}) + Jupiter ({jupiter_score}) + Saturn ({saturn_score}) + Mercurius ({mercurius_score}) + Smith ({smith_score}) + Concordia ({concordia_score}) = {total_score}")
