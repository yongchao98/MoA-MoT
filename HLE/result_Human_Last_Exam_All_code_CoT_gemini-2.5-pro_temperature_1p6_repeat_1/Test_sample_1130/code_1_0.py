import math

# --- Game State Information ---
# Player assets
houses_total = 12
houses_brick = 3
houses_tools = 3
provinces_with_houses = 7
colonists = 5
sestertii = 13
has_concordia_card = True

# Player cards
cards_vesta = 1
cards_jupiter = 2
cards_saturn = 2
cards_mercurius = 2
has_smith_card = True # A Minerva-type card
# The prompt does not mention any Mars cards, so we assume 0.
cards_mars = 0 

# Other game information
goods_produced_types = 4 # brick, food, tools, cloth

# --- Scoring Calculation ---
# 1. Vesta Score (1 VP per 10 sestertii, per card)
vesta_score = cards_vesta * math.floor(sestertii / 10)

# 2. Jupiter Score (1 VP per non-brick city, per card)
non_brick_houses = houses_total - houses_brick
jupiter_score = cards_jupiter * non_brick_houses

# 3. Saturn Score (1 VP per province with a house, per card)
saturn_score = cards_saturn * provinces_with_houses

# 4. Mercurius Score (2 VP per type of good produced, per card)
mercurius_score = cards_mercurius * goods_produced_types

# 5. Minerva Score (from Smith card: 2 VP per tool city)
# Smith card itself dictates the scoring, not a multiplier from another card.
if has_smith_card:
    minerva_score = 2 * houses_tools
else:
    minerva_score = 0

# 6. Mars Score (2 VP per colonist, per card)
mars_score = cards_mars * 2 * colonists

# 7. Concordia Card Score
if has_concordia_card:
    concordia_score = 7
else:
    concordia_score = 0

# --- Total Score ---
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + minerva_score + mars_score + concordia_score

# --- Final Output ---
print("The player's final score is calculated as follows:")
print(
    f"{vesta_score} (Vesta) + "
    f"{jupiter_score} (Jupiter) + "
    f"{saturn_score} (Saturn) + "
    f"{mercurius_score} (Mercurius) + "
    f"{minerva_score} (Minerva - Smith) + "
    f"{mars_score} (Mars) + "
    f"{concordia_score} (Concordia) = {total_score}"
)
print(f"\nFinal Score: {total_score}")
<<<54>>>