import math

# Game state variables
sestertii = 13
total_houses = 12
brick_houses = 3
provinces_with_houses = 7
goods_produced_types = 4 # brick, food, tool, cloth
tool_houses = 3
has_concordia_card = True

# Player's cards
num_vesta_cards = 1
num_jupiter_cards = 2
num_saturn_cards = 2
num_mercurius_cards = 2
num_smith_cards = 1 # This is a Minerva card for tools

# --- Score Calculation ---

# 1. Vesta Score
# 1 point for every 10 sestertii
vesta_score = num_vesta_cards * math.floor(sestertii / 10)
print(f"Vesta Score: {num_vesta_cards} card * floor({sestertii} sestertii / 10) = {vesta_score} points")

# 2. Jupiter Score
# 1 point for each non-brick city per card
non_brick_houses = total_houses - brick_houses
jupiter_score = num_jupiter_cards * non_brick_houses
print(f"Jupiter Score: {num_jupiter_cards} cards * ({total_houses} total houses - {brick_houses} brick houses) = {jupiter_score} points")

# 3. Saturn Score
# 2 points for each province with a house per card
saturn_score = num_saturn_cards * (provinces_with_houses * 2)
print(f"Saturn Score: {num_saturn_cards} cards * ({provinces_with_houses} provinces * 2) = {saturn_score} points")

# 4. Mercurius Score
# 2 points for each type of good produced per card
mercurius_score = num_mercurius_cards * (goods_produced_types * 2)
print(f"Mercurius Score: {num_mercurius_cards} cards * ({goods_produced_types} types of goods * 2) = {mercurius_score} points")

# 5. Minerva (Smith) Score
# 3 points for each tool-producing house
minerva_smith_score = num_smith_cards * (tool_houses * 3)
print(f"Minerva (Smith) Score: {num_smith_cards} card * ({tool_houses} tool houses * 3) = {minerva_smith_score} points")

# 6. Concordia Card Score
concordia_score = 7 if has_concordia_card else 0
print(f"Concordia Card Score: {concordia_score} points")

# 7. Total Score
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + minerva_smith_score + concordia_score
print("\n--- Final Score Calculation ---")
print(f"Total Score = {vesta_score} (Vesta) + {jupiter_score} (Jupiter) + {saturn_score} (Saturn) + {mercurius_score} (Mercurius) + {minerva_smith_score} (Smith) + {concordia_score} (Concordia)")
print(f"Final Score = {total_score}")

print(f'<<<{total_score}>>>')