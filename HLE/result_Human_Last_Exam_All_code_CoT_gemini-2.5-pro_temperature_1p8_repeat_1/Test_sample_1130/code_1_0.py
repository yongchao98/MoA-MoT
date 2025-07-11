def calculate_concordia_score():
    """
    Calculates the final score in a game of Concordia based on provided details.
    """
    # Given values
    sestertii = 13
    vesta_cards = 1
    
    jupiter_cards = 2
    total_provinces_with_houses = 7
    # Assumption: All 3 brick-producing houses are in a single province to maximize Jupiter score.
    # Therefore, the number of non-brick provinces with houses is total_provinces - 1.
    non_brick_provinces = total_provinces_with_houses - 1
    
    saturn_cards = 2
    
    mercurius_cards = 2
    good_types_produced = 4  # brick, food, tools, cloth

    smith_cards = 1
    tool_cities = 3
    smith_vp_per_city = 3

    concordia_card_vp = 7

    # Calculate points for each category
    vesta_points = (sestertii // 10) * vesta_cards
    jupiter_points = jupiter_cards * non_brick_provinces
    saturn_points = saturn_cards * total_provinces_with_houses
    mercurius_points = mercurius_cards * good_types_produced * 2
    minerva_points = smith_cards * tool_cities * smith_vp_per_city
    
    total_score = vesta_points + jupiter_points + saturn_points + mercurius_points + minerva_points + concordia_card_vp
    
    # Print the detailed scoring breakdown
    print("Final Score Calculation:")
    print(f"{vesta_points} (Vesta) + "
          f"{jupiter_points} (Jupiter) + "
          f"{saturn_points} (Saturnus) + "
          f"{mercurius_points} (Mercurius) + "
          f"{minerva_points} (Minerva - Smith) + "
          f"{concordia_card_vp} (Concordia) = {total_score}")

calculate_concordia_score()