# Player's assets
sestertii = 13
total_houses = 12
brick_houses = 3
tool_houses = 3
provinces = 7
# Number of different goods produced (brick, food, tools, cloth)
good_types_produced = 4

# Player's cards
num_vesta_cards = 1
num_jupiter_cards = 2
num_saturn_cards = 2
num_mercurius_cards = 2
num_smith_cards = 1
has_concordia_card = True

# Standard VP multipliers from Concordia rules
vesta_vp_multiplier = 1 # 1 VP per 10 sestertii
jupiter_vp_multiplier = 2 # 2 VP per non-brick city
saturn_vp_multiplier = 1 # 1 VP per province
mercurius_vp_multiplier = 2 # 2 VP per good type
smith_vp_multiplier = 3 # 3 VP per tool city
concordia_card_vp = 7

# --- Score Calculation ---

# 1. Vesta Score (Sestertii)
vesta_score = (sestertii // 10) * vesta_vp_multiplier * num_vesta_cards

# 2. Jupiter Score (Non-Brick Cities)
non_brick_cities = total_houses - brick_houses
jupiter_score = non_brick_cities * jupiter_vp_multiplier * num_jupiter_cards

# 3. Saturn Score (Provinces)
saturn_score = provinces * saturn_vp_multiplier * num_saturn_cards

# 4. Mercurius Score (Goods)
mercurius_score = good_types_produced * mercurius_vp_multiplier * num_mercurius_cards

# 5. Minerva Score (Smith - Tool Cities)
minerva_smith_score = tool_houses * smith_vp_multiplier * num_smith_cards

# 6. Concordia Card Score
concordia_score = concordia_card_vp if has_concordia_card else 0

# 7. Total Score
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + minerva_smith_score + concordia_score

# --- Output the results ---
print("Calculating Concordia Score:\n")
print(f"Vesta Points: ({sestertii} sestertii // 10) * {num_vesta_cards} Vesta Card(s) = {vesta_score}")
print(f"Jupiter Points: ({total_houses} total houses - {brick_houses} brick houses) * {jupiter_vp_multiplier} VP * {num_jupiter_cards} Jupiter Card(s) = {jupiter_score}")
print(f"Saturn Points: {provinces} provinces * {saturn_vp_multiplier} VP * {num_saturn_cards} Saturn Card(s) = {saturn_score}")
print(f"Mercurius Points: {good_types_produced} good types * {mercurius_vp_multiplier} VP * {num_mercurius_cards} Mercurius Card(s) = {mercurius_score}")
print(f"Minerva (Smith) Points: {tool_houses} tool houses * {smith_vp_multiplier} VP * {num_smith_cards} Smith Card(s) = {minerva_smith_score}")
print(f"Concordia Card Points: {concordia_score}\n")

print("--- Total Score ---")
print(f"{vesta_score} (Vesta) + {jupiter_score} (Jupiter) + {saturn_score} (Saturn) + {mercurius_score} (Mercurius) + {minerva_smith_score} (Minerva) + {concordia_score} (Concordia) = {total_score}")

print(f"\nThe player's total score is: {total_score}")
