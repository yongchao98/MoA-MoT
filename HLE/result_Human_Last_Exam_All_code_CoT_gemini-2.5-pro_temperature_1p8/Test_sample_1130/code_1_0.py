# Game assets
sestertii = 13
total_houses = 12
brick_houses = 3
provinces_with_houses = 7
goods_produced_types = 4  # Brick, Food, Tool, Cloth
colonists = 5

# Cards
vesta_cards = 1
jupiter_cards = 2
saturn_cards = 2
mercurius_cards = 2
mars_cards = 1  # Smith is a Mars card
has_concordia_card = True

# --- Score Calculation ---

# Vesta scoring: 1 point per 10 sestertii
vesta_score = (sestertii // 10) * vesta_cards
print(f"Vesta Points: {vesta_cards} card * ({sestertii} sestertii // 10) = {vesta_score}")

# Jupiter scoring: 1 point per non-brick city per card
non_brick_cities = total_houses - brick_houses
jupiter_score = non_brick_cities * jupiter_cards
print(f"Jupiter Points: {jupiter_cards} cards * {non_brick_cities} non-brick cities = {jupiter_score}")

# Saturn scoring: 1 point per province with at least one house per card
saturn_score = provinces_with_houses * saturn_cards
print(f"Saturn Points: {saturn_cards} cards * {provinces_with_houses} provinces = {saturn_score}")

# Mercurius scoring: 2 points per type of good produced per card
mercurius_score = mercurius_cards * goods_produced_types * 2
print(f"Mercurius Points: {mercurius_cards} cards * {goods_produced_types} types of goods * 2 = {mercurius_score}")

# Mars (Smith) scoring: 2 points per colonist per card
mars_score = colonists * 2 * mars_cards
print(f"Mars (Smith) Points: {mars_cards} card * {colonists} colonists * 2 = {mars_score}")

# Concordia card scoring: 7 points
concordia_score = 7 if has_concordia_card else 0
print(f"Concordia Card Points: {concordia_score}")

# Total score
total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + mars_score + concordia_score

# Final equation
print("\n--- Final Score Calculation ---")
print(f"{vesta_score} (Vesta) + {jupiter_score} (Jupiter) + {saturn_score} (Saturn) + {mercurius_score} (Mercurius) + {mars_score} (Mars) + {concordia_score} (Concordia) = {total_score}")
print(f"\nTotal Score: {total_score}")