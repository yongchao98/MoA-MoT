import math

# --- Player's Assets ---
# Cards
vesta_cards = 1
jupiter_cards = 2
saturn_cards = 2
mercurius_cards = 2
smith_cards = 1  # Smith is a Minerva card

# Game State
sestertii = 13
provinces_with_houses = 7
concordia_card_bonus = 7

# Production (houses in cities)
brick_houses = 3
food_houses = 4
tool_houses = 3
cloth_houses = 2

# --- Score Calculation ---

# 1. Vesta Score: 1 VP for every 10 sestertii per Vesta card
vesta_score = vesta_cards * math.floor(sestertii / 10)

# 2. Jupiter Score: 1 VP for each non-brick city per Jupiter card
non_brick_houses = food_houses + tool_houses + cloth_houses
jupiter_score = jupiter_cards * non_brick_houses

# 3. Saturn Score: 1 VP for each province with a house per Saturn card
saturn_score = saturn_cards * provinces_with_houses

# 4. Mercurius Score: 2 VP for each type of good produced per Mercurius card
# The player produces 4 types of goods: brick, food, tools, cloth
produced_goods_types = 4
mercurius_score = mercurius_cards * produced_goods_types * 2 # Each Mercurius card grants 2 VP per good type

# 5. Smith (Minerva) Score: 2 VP for each tool-producing city
smith_score = tool_houses * 2

# 6. Concordia Card Score
concordia_score = concordia_card_bonus

# --- Total Score ---
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + smith_score + concordia_score

# --- Final Output ---
print("The final score is calculated as follows:")
print("Vesta + Jupiter + Saturn + Mercurius + Smith + Concordia Card = Total")
print(f"{vesta_score} + {jupiter_score} + {saturn_score} + {mercurius_score} + {smith_score} + {concordia_score} = {total_score}")

print("The score breakdown is:")
print(f"- Vesta ({vesta_cards} card * floor({sestertii}/10) sestertii): {vesta_score} points")
print(f"- Jupiter ({jupiter_cards} cards * {non_brick_houses} non-brick cities): {jupiter_score} points")
print(f"- Saturn ({saturn_cards} cards * {provinces_with_houses} provinces): {saturn_score} points")
print(f"- Mercurius ({mercurius_cards} cards * {produced_goods_types} good types * 2): {mercurius_score} points")
print(f"- Smith (Minerva) ({tool_houses} tool cities * 2): {smith_score} points")
print(f"- Concordia Card: {concordia_score} points")
<<<73>>>