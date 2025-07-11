def calculate_aog_score():
    """
    Calculates the final score for a game of Age of Galaxy based on provided details.
    """
    # Player's assets
    credits = 7
    prestige = 5
    relics = 2
    cruisers = 5
    colonized_planets = 3
    developed_planets = 1
    conquered_planets = 1
    # Only green (Exploration) techs count for the Us'ud bonus
    num_green_techs = 2 

    # --- VP Calculations ---

    # 1. Base VP
    prestige_vp = prestige  # 1 VP per Prestige
    relic_vp = relics * 3  # 3 VP per Relic
    colonized_vp_base = colonized_planets * 1 # 1 VP per Colonized Planet
    developed_vp_base = developed_planets * 3 # 3 VP per Developed Planet
    conquered_vp = conquered_planets * 2  # 2 VP per Conquered Planet

    # 2. Alliance VP
    # Chaeilki: 1 VP per 2 credits
    chaeilki_vp = credits // 2
    # Humans: 2 VP per colonized planet
    humans_vp = colonized_planets * 2
    # Us'ud: 2 VP per green tech
    usud_vp = num_green_techs * 2

    # 3. Ideology VP (Legarchaea)
    # 2 VP per developed planet
    legarchaea_dev_vp = developed_planets * 2
    # 1 VP per 2 cruisers
    legarchaea_cruiser_vp = cruisers // 2

    # 4. Sum Total
    total_vp = (prestige_vp + relic_vp + colonized_vp_base + developed_vp_base +
                conquered_vp + chaeilki_vp + humans_vp + usud_vp +
                legarchaea_dev_vp + legarchaea_cruiser_vp)
    
    # --- Output ---
    print("VP Breakdown and Calculation:\n")

    # The prompt asks for each number in the final equation to be output.
    # The following print statements show the origin of each number.
    print(f"Prestige VP: {prestige_vp}")
    print(f"Relic VP: {relics} * 3 = {relic_vp}")
    print(f"Base Colonized Planet VP: {colonized_planets} * 1 = {colonized_vp_base}")
    print(f"Base Developed Planet VP: {developed_planets} * 3 = {developed_vp_base}")
    print(f"Conquered Planet VP: {conquered_planets} * 2 = {conquered_vp}")
    print(f"Chaeilki Alliance VP (from credits): {credits} // 2 = {chaeilki_vp}")
    print(f"Humans Alliance VP (from colonized planets): {colonized_planets} * 2 = {humans_vp}")
    print(f"Us'ud Alliance VP (from green techs): {num_green_techs} * 2 = {usud_vp}")
    print(f"Legarchaea Ideology VP (from developed planets): {developed_planets} * 2 = {legarchaea_dev_vp}")
    print(f"Legarchaea Ideology VP (from cruisers): {cruisers} // 2 = {legarchaea_cruiser_vp}")
    
    print("-" * 35)
    
    final_equation = (f"{prestige_vp} + {relic_vp} + {colonized_vp_base} + "
                      f"{developed_vp_base} + {conquered_vp} + {chaeilki_vp} + "
                      f"{humans_vp} + {usud_vp} + {legarchaea_dev_vp} + "
                      f"{legarchaea_cruiser_vp}")

    print(f"Final VP Calculation:\n{final_equation} = {total_vp}")


calculate_aog_score()
<<<36>>>