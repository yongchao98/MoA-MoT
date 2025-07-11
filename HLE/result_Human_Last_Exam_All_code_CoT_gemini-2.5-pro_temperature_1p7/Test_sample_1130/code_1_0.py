import math

# Step 1: Define the player's game state from the problem description.
# Assets
provinces_with_houses = 7
sestertii = 13
has_concordia_card = True

# Houses by production type
total_houses = 12
brick_houses = 3
food_houses = 4
tool_houses = 3
cloth_houses = 2
wine_houses = 0  # Not mentioned, assumed 0

# Goods in storehouse
storehouse_cloth = 1
storehouse_tools = 4
storehouse_brick = 1

# Personality Cards
vesta_cards = 1
saturn_cards = 2
jupiter_cards = 2
mercurius_cards = 2
smith_cards = 1  # Smith is a specialist card for Tools

# Step 2: Calculate the score for each category.

# Vesta Score Calculation
# Assumption: Since goods values are not given, standard values are used: Brick=3, Tools=5, Cloth=6.
goods_value = {"brick": 3, "tools": 5, "cloth": 6}
goods_total_value = (storehouse_brick * goods_value["brick"]) + \
                    (storehouse_tools * goods_value["tools"]) + \
                    (storehouse_cloth * goods_value["cloth"])
total_monetary_value = sestertii + goods_total_value
vesta_score = math.floor(total_monetary_value / 10)

# Saturn Score Calculation (1 VP per province per card)
saturn_score = saturn_cards * provinces_with_houses

# Jupiter Score Calculation (1 VP per non-brick house per card)
non_brick_houses = total_houses - brick_houses
jupiter_score = jupiter_cards * non_brick_houses

# Mercurius Score Calculation (2 VPs per type of good produced per card)
goods_types_produced = 0
if brick_houses > 0: goods_types_produced += 1
if food_houses > 0: goods_types_produced += 1
if tool_houses > 0: goods_types_produced += 1
if cloth_houses > 0: goods_types_produced += 1
if wine_houses > 0: goods_types_produced += 1
mercurius_score = mercurius_cards * goods_types_produced * 2

# Specialist Card (Smith) Score Calculation (4 VPs per Tool house)
smith_score_per_house = 4
smith_score = tool_houses * smith_score_per_house

# Concordia Card Score
concordia_score = 7 if has_concordia_card else 0

# Total Score Calculation
total_score = vesta_score + saturn_score + jupiter_score + mercurius_score + smith_score + concordia_score

# Step 3: Print the results in a clear, step-by-step format.
print("Concordia Score Calculation Breakdown:")
print("-------------------------------------")
print(f"1. Vesta Score (based on {vesta_cards} card):")
print(f"   Value of goods ({storehouse_brick} brick, {storehouse_tools} tools, {storehouse_cloth} cloth) is {goods_total_value}.")
print(f"   Total value ({sestertii} sestertii + {goods_total_value} goods) = {total_monetary_value}.")
print(f"   Score = floor({total_monetary_value} / 10) = {vesta_score} points.")

print(f"\n2. Saturn Score (based on {saturn_cards} cards):")
print(f"   {saturn_cards} cards * {provinces_with_houses} provinces = {saturn_score} points.")

print(f"\n3. Jupiter Score (based on {jupiter_cards} cards):")
print(f"   {jupiter_cards} cards * {non_brick_houses} non-brick houses = {jupiter_score} points.")

print(f"\n4. Mercurius Score (based on {mercurius_cards} cards):")
print(f"   {mercurius_cards} cards * {goods_types_produced} types of goods * 2 VP/type = {mercurius_score} points.")

print(f"\n5. Smith Score (based on {smith_cards} card):")
print(f"   {tool_houses} tool houses * {smith_score_per_house} VP/house = {smith_score} points.")

print(f"\n6. Concordia Card Score:")
print(f"   Bonus for owning the card = {concordia_score} points.")

print("\n-------------------------------------")
print("Final Score Equation:")
print(f"{vesta_score} (Vesta) + {saturn_score} (Saturn) + {jupiter_score} (Jupiter) + {mercurius_score} (Mercurius) + {smith_score} (Smith) + {concordia_score} (Concordia) = {total_score}")
print("-------------------------------------")
print(f"The player's total score is: {total_score}")