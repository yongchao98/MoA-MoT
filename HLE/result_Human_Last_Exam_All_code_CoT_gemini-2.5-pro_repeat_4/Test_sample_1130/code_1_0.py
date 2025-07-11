# Player's assets and cards as described
houses_brick = 3
houses_food = 4
houses_tools = 3
houses_cloth = 2
total_houses = 12
provinces = 7
sestertii = 13
has_concordia_card = True

cards_saturn = 2
cards_jupiter = 2
cards_vesta = 1
cards_mercurius = 2
# The Smith card is a Minerva card for Tools.

# --- Score Calculation Step-by-Step ---

# Vesta Score: 1 VP for every 10 sestertii
vesta_score = (sestertii // 10) * cards_vesta

# Jupiter Score: 1 VP per non-brick city, per card
non_brick_cities = total_houses - houses_brick
jupiter_score = cards_jupiter * non_brick_cities

# Saturn Score: 1 VP per province with a house, per card
saturn_score = cards_saturn * provinces

# Mercurius Score: 2 VP for each type of good produced, per card
# The player produces 4 types of goods: bricks, food, tools, and cloth.
num_good_types = 4
mercurius_score = cards_mercurius * (num_good_types * 2)

# Minerva (Smith) Score: 3 VP for each tool-producing city
minerva_smith_score = houses_tools * 3

# Concordia Card Score: 7 VP for ending the game
concordia_score = 7 if has_concordia_card else 0

# --- Total Score ---
all_scores = [
    vesta_score,
    jupiter_score,
    saturn_score,
    mercurius_score,
    minerva_smith_score,
    concordia_score
]
total_score = sum(all_scores)

# --- Output the results ---
print("Score Breakdown:")
print(f"Vesta: {vesta_score} points")
print(f"Jupiter: {jupiter_score} points")
print(f"Saturn: {saturn_score} points")
print(f"Mercurius: {mercurius_score} points")
print(f"Minerva (Smith): {minerva_smith_score} points")
print(f"Concordia Card: {concordia_score} points")

# Print the final equation with each individual score value
final_equation_str = " + ".join(map(str, all_scores))
print("\nFinal Score Calculation:")
print(f"{final_equation_str} = {total_score}")