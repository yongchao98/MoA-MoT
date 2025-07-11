# Game state variables
# Houses by production type
brick_houses = 3
food_houses = 4
tool_houses = 3
cloth_houses = 2
total_houses = brick_houses + food_houses + tool_houses + cloth_houses

# Other assets
provinces = 7
colonists = 5
sestertii = 13
has_concordia_card = True

# Storehouse contents
storehouse_cloth = 1
storehouse_tools = 4
storehouse_brick = 1

# Player's cards
vesta_cards = 1
jupiter_cards = 2
saturn_cards = 2
mercurius_cards = 2
smith_cards = 1

# --- Score Calculation ---

# 1. Vesta Score
good_values = {'brick': 3, 'food': 4, 'tool': 5, 'cloth': 6}
storehouse_value = (storehouse_brick * good_values['brick']) + \
                   (storehouse_tools * good_values['tool']) + \
                   (storehouse_cloth * good_values['cloth'])
total_assets_value = sestertii + storehouse_value
vesta_points = (total_assets_value // 10) * vesta_cards

# 2. Jupiter Score
non_brick_houses = total_houses - brick_houses
jupiter_points = non_brick_houses * jupiter_cards

# 3. Saturn Score
saturn_points = provinces * saturn_cards

# 4. Mercurius Score
# The player produces bricks, food, tools, and cloth.
produced_good_types = 4
mercurius_points = (produced_good_types * 2) * mercurius_cards

# 5. Minerva (Smith) Score
# Smith card gives 4 points per tool-producing city
smith_points = tool_houses * 4 * smith_cards

# 6. Concordia Card Score
concordia_points = 7 if has_concordia_card else 0

# 7. Total Score
total_score = vesta_points + jupiter_points + saturn_points + mercurius_points + smith_points + concordia_points

# --- Output the final equation ---
print(f"The final score is calculated as follows:")
print(f"Vesta ({vesta_cards} card): {vesta_points}")
print(f"Jupiter ({jupiter_cards} cards): {jupiter_points}")
print(f"Saturn ({saturn_cards} cards): {saturn_points}")
print(f"Mercurius ({mercurius_cards} cards): {mercurius_points}")
print(f"Smith ({smith_cards} card): {smith_points}")
print(f"Concordia Card: {concordia_points}")
print("-" * 20)
print(f"Total Score = {vesta_points} + {jupiter_points} + {saturn_points} + {mercurius_points} + {smith_points} + {concordia_points} = {total_score}")
print(f"<<<{total_score}>>>")