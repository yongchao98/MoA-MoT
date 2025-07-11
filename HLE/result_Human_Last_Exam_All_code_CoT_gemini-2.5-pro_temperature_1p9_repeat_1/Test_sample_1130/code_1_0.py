import math

def calculate_concordia_score():
    """
    Calculates and prints the Concordia game score based on provided game state.
    """
    # --- Player's Assets ---
    num_brick_houses = 3
    num_food_houses = 4
    num_tool_houses = 3
    num_cloth_houses = 2
    total_houses = num_brick_houses + num_food_houses + num_tool_houses + num_cloth_houses
    
    num_provinces_with_houses = 7
    num_colonists = 5
    sestertii = 13
    has_concordia_card = True
    
    # --- Player's Cards ---
    num_vesta_cards = 1
    num_jupiter_cards = 2
    num_saturn_cards = 2
    num_mercurius_cards = 2
    num_smith_cards = 1  # Minerva card for Tools
    # Assuming the player has the default starting Mars card.
    num_mars_cards = 1
    
    # --- Score Calculation ---
    total_score = 0
    final_equation_parts = []
    
    print("--- Concordia Score Breakdown ---")

    # 1. Vesta Score (Sestertii)
    vesta_points_per_card = sestertii // 10
    vesta_score = vesta_points_per_card * num_vesta_cards
    total_score += vesta_score
    final_equation_parts.append(str(vesta_score))
    print(f"Vesta Score (Sestertii): {num_vesta_cards} card * {vesta_points_per_card} points = {vesta_score} VPs")
    
    # 2. Jupiter Score (Non-Brick Cities)
    num_non_brick_cities = total_houses - num_brick_houses
    jupiter_score = num_non_brick_cities * num_jupiter_cards
    total_score += jupiter_score
    final_equation_parts.append(str(jupiter_score))
    print(f"Jupiter Score (Cities): {num_jupiter_cards} cards * {num_non_brick_cities} non-brick cities = {jupiter_score} VPs")
    
    # 3. Saturn Score (Provinces)
    saturn_score = num_provinces_with_houses * num_saturn_cards
    total_score += saturn_score
    final_equation_parts.append(str(saturn_score))
    print(f"Saturn Score (Provinces): {num_saturn_cards} cards * {num_provinces_with_houses} provinces = {saturn_score} VPs")
    
    # 4. Mercurius Score (Goods Production)
    num_good_types_produced = 4  # Brick, Food, Tools, Cloth
    points_per_good_type = 2
    mercurius_score_per_card = num_good_types_produced * points_per_good_type
    mercurius_score = mercurius_score_per_card * num_mercurius_cards
    total_score += mercurius_score
    final_equation_parts.append(str(mercurius_score))
    print(f"Mercurius Score (Goods): {num_mercurius_cards} cards * ({num_good_types_produced} types * {points_per_good_type}) = {mercurius_score} VPs")
    
    # 5. Minerva Score (Smith card for Tools)
    points_per_tool_city = 3
    minerva_score = num_tool_houses * points_per_tool_city * num_smith_cards
    total_score += minerva_score
    final_equation_parts.append(str(minerva_score))
    print(f"Minerva Score (Smith): {num_smith_cards} card * ({num_tool_houses} tool cities * {points_per_tool_city}) = {minerva_score} VPs")
    
    # 6. Mars Score (Colonists)
    points_per_colonist = 2
    mars_score = num_colonists * points_per_colonist * num_mars_cards
    total_score += mars_score
    final_equation_parts.append(str(mars_score))
    print(f"Mars Score (Colonists): {num_mars_cards} card * ({num_colonists} colonists * {points_per_colonist}) = {mars_score} VPs")
    
    # 7. Concordia Card
    concordia_score = 7 if has_concordia_card else 0
    total_score += concordia_score
    final_equation_parts.append(str(concordia_score))
    print(f"Concordia Card Bonus: {concordia_score} VPs")
    
    # --- Final Score ---
    print("\n--- Final Score ---")
    final_equation_str = f"{final_equation_parts[0]} (Vesta) + {final_equation_parts[1]} (Jupiter) + {final_equation_parts[2]} (Saturn) + {final_equation_parts[3]} (Mercurius) + {final_equation_parts[4]} (Minerva) + {final_equation_parts[5]} (Mars) + {final_equation_parts[6]} (Concordia)"
    print(f"Total Score = {final_equation_str} = {total_score} VPs")
    return total_score

if __name__ == '__main__':
    final_score = calculate_concordia_score()
    print(f"\n<<<The final score is: {final_score}>>>")
