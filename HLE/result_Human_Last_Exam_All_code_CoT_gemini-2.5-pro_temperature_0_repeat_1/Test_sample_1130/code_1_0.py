# Game state variables
sestertii = 13
total_houses = 12
brick_houses = 3
provinces_with_houses = 7
# Goods produced: brick, food, tool, cloth
num_good_types = 4
tool_houses = 3
has_concordia_card = True

# Card counts
vesta_cards = 1
jupiter_cards = 2
saturn_cards = 2
mercurius_cards = 2
smith_cards = 1 # This is a Minerva card

# --- Score Calculation ---

# Vesta score (1 point per 10 sestertii)
vesta_score = (sestertii // 10) * vesta_cards

# Jupiter score (1 point per non-brick city)
non_brick_houses = total_houses - brick_houses
jupiter_score = non_brick_houses * jupiter_cards

# Saturn score (1 point per province with a house)
saturn_score = provinces_with_houses * saturn_cards

# Mercurius score (2 points per type of good produced)
mercurius_score = num_good_types * 2 * mercurius_cards

# Minerva (Smith) score (3 points per tool city)
smith_score = tool_houses * 3 * smith_cards

# Concordia card score
concordia_score = 7 if has_concordia_card else 0

# Total score
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + smith_score + concordia_score

# --- Output ---
print("Calculating the final score:")
print(f"Vesta Score: {vesta_score}")
print(f"Jupiter Score: {jupiter_score}")
print(f"Saturn Score: {saturn_score}")
print(f"Mercurius Score: {mercurius_score}")
print(f"Smith (Minerva) Score: {smith_score}")
print(f"Concordia Card Score: {concordia_score}")
print("\nFinal Score Calculation:")
print(f"{vesta_score} (Vesta) + {jupiter_score} (Jupiter) + {saturn_score} (Saturn) + {mercurius_score} (Mercurius) + {smith_score} (Smith) + {concordia_score} (Concordia) = {total_score}")
