def calculate_concordia_score():
    """
    Calculates the final score in a game of Concordia based on the provided state.
    """
    # --- Game State Initialization ---
    sestertii = 13
    total_houses = 12
    provinces_with_houses = 7

    # Houses by production type
    brick_houses = 3
    food_houses = 4
    tool_houses = 3
    cloth_houses = 2
    
    # Player's cards
    vesta_cards = 1
    jupiter_cards = 2
    saturn_cards = 2
    mercurius_cards = 2
    smith_card = 1 # This is a Minerva card for Tools
    has_concordia_card = True
    
    # --- Scoring Calculation ---

    # 1. Vesta Score (based on Sestertii)
    # 1 VP for every 10 sestertii.
    vesta_score = sestertii // 10
    
    # 2. Jupiter Score (non-brick cities)
    # VP per card = total houses - brick houses
    non_brick_houses = total_houses - brick_houses
    jupiter_score = jupiter_cards * non_brick_houses
    
    # 3. Saturn Score (provinces)
    # VP per card = number of provinces with at least one house
    saturn_score = saturn_cards * provinces_with_houses
    
    # 4. Mercurius Score (types of goods produced)
    # 2 VP per type of good produced, per card
    goods_types_produced = 0
    if brick_houses > 0:
        goods_types_produced += 1
    if food_houses > 0:
        goods_types_produced += 1
    if tool_houses > 0:
        goods_types_produced += 1
    if cloth_houses > 0:
        goods_types_produced += 1
    mercurius_score = mercurius_cards * goods_types_produced * 2

    # 5. Minerva Score (Smith card for Tool cities)
    # Smith is a Minerva card for Tools. Tools are worth 5 VP per city.
    tool_city_value = 5
    minerva_score = smith_card * tool_houses * tool_city_value
    
    # 6. Concordia Card Score
    concordia_score = 7 if has_concordia_card else 0
    
    # --- Total Score ---
    total_score = vesta_score + jupiter_score + saturn_score + mercurius_score + minerva_score + concordia_score

    # --- Print Results ---
    print("Concordia Score Calculation:")
    print(f"- Vesta Score (from 13 sestertii): {vesta_score} VP")
    print(f"- Jupiter Score (2 cards * {non_brick_houses} non-brick cities): {jupiter_score} VP")
    print(f"- Saturn Score (2 cards * {provinces_with_houses} provinces): {saturn_score} VP")
    print(f"- Mercurius Score (2 cards * {goods_types_produced} goods types * 2): {mercurius_score} VP")
    print(f"- Minerva (Smith) Score ({tool_houses} tool cities * {tool_city_value} value): {minerva_score} VP")
    print(f"- Concordia Card Score: {concordia_score} VP")
    print("-" * 20)
    print("Final Score Equation:")
    print(f"Total Score = {vesta_score} + {jupiter_score} + {saturn_score} + {mercurius_score} + {minerva_score} + {concordia_score} = {total_score}")

calculate_concordia_score()
<<<71>>>