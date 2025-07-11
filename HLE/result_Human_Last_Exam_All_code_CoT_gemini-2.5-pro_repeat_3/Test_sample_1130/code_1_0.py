# Game state variables based on the user's description
sestertii = 13
provinces_with_houses = 7
colonists = 5
food_houses = 4
tool_houses = 3
cloth_houses = 2

# Card counts from the player's hand
vesta_cards = 1
jupiter_cards = 2
saturn_cards = 2
mercurius_cards = 2
smith_cards = 1  # A Minerva card for Tools
mars_cards = 0   # No Mars cards were listed

# Fixed bonus for the Concordia card
concordia_card_bonus = 7

# --- Score Calculation ---

# 1. Vesta Score (for sestertii)
# Formula: num_cards * floor(sestertii / 10)
vesta_value = sestertii // 10
vesta_score = vesta_cards * vesta_value

# 2. Jupiter Score (for non-brick cities)
# Formula: num_cards * num_non_brick_cities
non_brick_cities = food_houses + tool_houses + cloth_houses
jupiter_score = jupiter_cards * non_brick_cities

# 3. Saturn Score (for provinces)
# Formula: num_cards * num_provinces
saturn_score = saturn_cards * provinces_with_houses

# 4. Mercurius Score (for types of goods)
# Goods produced are brick, food, tools, and cloth (4 types)
goods_types_produced = 4
mercurius_score = mercurius_cards * goods_types_produced

# 5. Minerva Score (for Smith - tool cities)
# Formula: num_smith_cards * num_tool_cities
minerva_smith_score = smith_cards * tool_houses

# 6. Mars Score (for colonists)
# Formula: num_cards * num_colonists
# The player has no Mars cards, so the score is 0.
mars_score = mars_cards * colonists

# 7. Concordia Card Score
concordia_score = concordia_card_bonus

# Calculate the total score
total_score = (vesta_score + jupiter_score + saturn_score +
               mercurius_score + minerva_smith_score + mars_score +
               concordia_score)

# --- Output the results ---
print("Concordia Score Breakdown:")
print(f"Vesta:       {vesta_cards} card * floor({sestertii}/10) = {vesta_score} points")
print(f"Jupiter:     {jupiter_cards} cards * {non_brick_cities} cities = {jupiter_score} points")
print(f"Saturn:      {saturn_cards} cards * {provinces_with_houses} provinces = {saturn_score} points")
print(f"Mercurius:   {mercurius_cards} cards * {goods_types_produced} goods types = {mercurius_score} points")
print(f"Minerva:     {smith_cards} Smith card * {tool_houses} tool cities = {minerva_smith_score} points")
print(f"Mars:        {mars_cards} cards * {colonists} colonists = {mars_score} points")
print(f"Concordia:   Bonus for Concordia Card = {concordia_score} points")
print("-" * 30)

# Print the final equation as requested
print("Final Score Equation:")
print(f"{vesta_score} + {jupiter_score} + {saturn_score} + {mercurius_score} + {minerva_smith_score} + {mars_score} + {concordia_score} = {total_score}")
<<<51>>>