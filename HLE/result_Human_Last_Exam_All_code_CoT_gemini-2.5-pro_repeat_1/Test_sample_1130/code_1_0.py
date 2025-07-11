# Game state variables
num_houses_brick = 3
num_houses_food = 4
num_houses_tool = 3
num_houses_cloth = 2
num_provinces = 7
sestertii = 13
has_concordia_card = True
storehouse_cloth = 1
storehouse_tools = 4
storehouse_brick = 1
num_saturn_cards = 2
num_jupiter_cards = 2
num_vesta_cards = 1
num_mercurius_cards = 2
has_smith_card = True

# Scoring logic

# 1. Vesta Score
vesta_goods_value = (storehouse_cloth * 7) + (storehouse_tools * 5) + (storehouse_brick * 3)
vesta_total_value = vesta_goods_value + sestertii
vesta_score = vesta_total_value // 10
print(f"Vesta Score: {vesta_score} points from (({storehouse_cloth}*7 + {storehouse_tools}*5 + {storehouse_brick}*3) + {sestertii}) // 10")

# 2. Jupiter Score
num_non_brick_cities = num_houses_food + num_houses_tool + num_houses_cloth
jupiter_score = num_jupiter_cards * num_non_brick_cities
print(f"Jupiter Score: {jupiter_score} points from {num_jupiter_cards} * ({num_houses_food} + {num_houses_tool} + {num_houses_cloth})")

# 3. Saturn Score
saturn_score = num_saturn_cards * num_provinces
print(f"Saturn Score: {saturn_score} points from {num_saturn_cards} * {num_provinces}")

# 4. Mercurius Score
num_good_types = 0
if num_houses_brick > 0: num_good_types += 1
if num_houses_food > 0: num_good_types += 1
if num_houses_tool > 0: num_good_types += 1
if num_houses_cloth > 0: num_good_types += 1
# In this case, we know it's 4 types.
mercurius_score = num_mercurius_cards * num_good_types * 2
print(f"Mercurius Score: {mercurius_score} points from {num_mercurius_cards} * {num_good_types} * 2")

# 5. Minerva (Smith) Score
minerva_score = 0
if has_smith_card:
    # Smith is for tools and scores 4 points per tool city
    smith_score = num_houses_tool * 4
    minerva_score += smith_score
    print(f"Minerva (Smith) Score: {minerva_score} points from {num_houses_tool} * 4")

# 6. Concordia Card Score
concordia_score = 7 if has_concordia_card else 0
print(f"Concordia Card Score: {concordia_score} points")

# 7. Total Score
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + minerva_score + concordia_score
print("\n--- Final Score Calculation ---")
print(f"Total Score = {vesta_score} (Vesta) + {jupiter_score} (Jupiter) + {saturn_score} (Saturn) + {mercurius_score} (Mercurius) + {minerva_score} (Smith) + {concordia_score} (Concordia)")
print(f"Final Score = {total_score}")

<<<71>>>