import math

# --- Game State Variables ---

# Houses and Production
brick_cities = 3
food_cities = 4
tool_cities = 3
cloth_cities = 2
total_houses = 12 # 3 + 4 + 3 + 2

# Other Assets
provinces = 7
sestertii = 13
has_concordia_card = True

# Player Cards
vesta_cards = 1
jupiter_cards = 2
saturn_cards = 2
mercurius_cards = 2
smith_cards = 1 # This is a Minerva card for Tools

# --- Score Calculation ---

# 1. Vesta Score (for sestertii)
# Each Vesta card gives 1 VP for every 10 sestertii.
vesta_score = (sestertii // 10) * vesta_cards
print(f"Vesta Score: For {vesta_cards} Vesta card, the score is ({sestertii} sestertii // 10) * {vesta_cards} = {vesta_score} points")

# 2. Jupiter Score (for non-brick cities)
# Each Jupiter card gives 1 VP for each non-brick city.
non_brick_cities = food_cities + tool_cities + cloth_cities
jupiter_score = non_brick_cities * jupiter_cards
print(f"Jupiter Score: For {jupiter_cards} Jupiter cards, the score is {non_brick_cities} non-brick cities * {jupiter_cards} = {jupiter_score} points")

# 3. Saturn Score (for provinces)
# Each Saturn card gives 1 VP for each province with at least one house.
saturn_score = provinces * saturn_cards
print(f"Saturn Score: For {saturn_cards} Saturn cards, the score is {provinces} provinces * {saturn_cards} = {saturn_score} points")

# 4. Mercurius Score (for types of goods produced)
# Each Mercurius card gives 2 VP for each type of good produced.
goods_produced = ["brick", "food", "tool", "cloth"]
num_good_types = len(goods_produced)
vp_per_good_type = 2
mercurius_score = num_good_types * vp_per_good_type * mercurius_cards
print(f"Mercurius Score: For {mercurius_cards} Mercurius cards, the score is {num_good_types} good types * {vp_per_good_type} VP * {mercurius_cards} = {mercurius_score} points")

# 5. Minerva (Smith) Score (for tool cities)
# The Smith card gives 5 VP for each city producing Tools.
vp_per_tool_city = 5
minerva_smith_score = tool_cities * vp_per_tool_city
print(f"Minerva (Smith) Score: The score is {tool_cities} tool cities * {vp_per_tool_city} VP = {minerva_smith_score} points")

# 6. Concordia Card Score
concordia_score = 7 if has_concordia_card else 0
print(f"Concordia Card Score: The Concordia card is worth {concordia_score} points")

# --- Total Score ---
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + minerva_smith_score + concordia_score

print("\n--- Final Score Calculation ---")
print(f"Total Score = {vesta_score} (Vesta) + {jupiter_score} (Jupiter) + {saturn_score} (Saturn) + {mercurius_score} (Mercurius) + {minerva_smith_score} (Smith) + {concordia_score} (Concordia)")
print(f"The player's total score is: {total_score}")

# The final answer in the requested format
# print(f"<<<{total_score}>>>")