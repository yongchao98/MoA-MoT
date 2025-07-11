# Game state variables
sestertii = 13
total_houses = 12
brick_houses = 3
tool_houses = 3
provinces = 7
goods_types_produced = 4 # Brick, Food, Tools, Cloth

# Card counts
vesta_cards = 1
jupiter_cards = 2
saturn_cards = 2
mercurius_cards = 2
# Specialist cards
has_smith_card = True
# Game end card
has_concordia_card = True

# --- Score Calculation ---

# Vesta Score (based on sestertii)
# 1 point per 10 sestertii, multiplied by the number of Vesta cards
vesta_score = (sestertii // 10) * vesta_cards

# Jupiter Score (based on non-brick cities)
non_brick_cities = total_houses - brick_houses
jupiter_score = non_brick_cities * jupiter_cards

# Saturn Score (based on provinces)
saturn_score = provinces * saturn_cards

# Mercurius Score (based on types of goods produced)
mercurius_score = goods_types_produced * mercurius_cards

# Smith Score (specialist card for tool cities)
smith_score = 0
if has_smith_card:
    smith_score = tool_houses * 1 # 1 VP per tool city

# Concordia Card Score
concordia_score = 0
if has_concordia_card:
    concordia_score = 7

# --- Total Score ---
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + smith_score + concordia_score

# --- Print the final equation ---
print("The player's final score is calculated as follows:")
print(f"{vesta_score} (from Vesta) + {jupiter_score} (from Jupiter) + {saturn_score} (from Saturn) + {mercurius_score} (from Mercurius) + {smith_score} (from Smith) + {concordia_score} (from Concordia Card) = {total_score}")

print("\nFinal Score:")
print(total_score)
<<<51>>>