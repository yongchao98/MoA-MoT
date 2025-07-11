def calculate_aog_vp():
    """
    Calculates the total Victory Points (VP) for a player in Age of Galaxy.
    """
    # 1. Alliance VP
    chaeilki_vp = 5
    humans_vp = 4
    usud_vp = 4
    alliance_bonus_vp = 3  # For having 3 races
    total_alliance_vp = chaeilki_vp + humans_vp + usud_vp + alliance_bonus_vp

    # 2. Ideology VP
    legarchaea_vp = 6

    # 3. Resource VP
    # VP from Credits, Productivity, Discovery, Influence (1 VP per 5)
    credits = 7
    productivity = 5
    discovery = 4
    influence = 3
    basic_resources_total = credits + productivity + discovery + influence
    basic_resources_vp = basic_resources_total // 5

    # VP from Prestige (1 VP per 1)
    prestige = 5
    prestige_vp = prestige

    # VP from Relics (3 VP per 1)
    relics = 2
    relics_vp = relics * 3

    # VP from Cruisers (1 VP per 1)
    cruisers = 5
    cruisers_vp = cruisers * 1
    
    total_resource_vp = basic_resources_vp + prestige_vp + relics_vp + cruisers_vp

    # 4. Planet VP
    # VP from Colonized Planets (1 VP per 1)
    colonized_planets = 3
    colonized_vp = colonized_planets * 1

    # VP from Developed Planets (3 VP per 1)
    developed_planets = 1
    developed_vp = developed_planets * 3

    # VP from Conquered Planets (2 VP per 1)
    conquered_planets = 1
    conquered_vp = conquered_planets * 2

    total_planet_vp = colonized_vp + developed_vp + conquered_vp
    
    # 5. Research VP (n*(n+1)/2 for n technologies)
    num_techs = 3
    research_vp = num_techs * (num_techs + 1) // 2
    
    # 6. Total VP
    total_vp = total_alliance_vp + legarchaea_vp + total_resource_vp + total_planet_vp + research_vp

    # Print the breakdown
    print("Calculating Total Victory Points (VP):")
    print("\n--- Alliance VP ---")
    print(f"Chaeilki: {chaeilki_vp} VP")
    print(f"Humans: {humans_vp} VP")
    print(f"Us'ud: {usud_vp} VP")
    print(f"3-Race Alliance Bonus: {alliance_bonus_vp} VP")
    
    print("\n--- Ideology VP ---")
    print(f"Legarchaea: {legarchaea_vp} VP")
    
    print("\n--- Resource VP ---")
    print(f"Basic Resources ({basic_resources_total} total): {basic_resources_vp} VP")
    print(f"Prestige ({prestige} total): {prestige_vp} VP")
    print(f"Relics ({relics} total): {relics_vp} VP")
    print(f"Cruisers ({cruisers} total): {cruisers_vp} VP")

    print("\n--- Planet VP ---")
    print(f"Colonized Planets ({colonized_planets} total): {colonized_vp} VP")
    print(f"Developed Planets ({developed_planets} total): {developed_vp} VP")
    print(f"Conquered Planets ({conquered_planets} total): {conquered_vp} VP")
    
    print("\n--- Research VP ---")
    print(f"Technologies ({num_techs} total): {research_vp} VP")

    print("\n--- Final Calculation ---")
    print(f"Total VP = {chaeilki_vp} (Chaeilki) + {humans_vp} (Humans) + {usud_vp} (Us'ud) + {alliance_bonus_vp} (Alliance) + {legarchaea_vp} (Ideology) + {basic_resources_vp} (Resources) + {prestige_vp} (Prestige) + {relics_vp} (Relics) + {cruisers_vp} (Cruisers) + {colonized_vp} (Colonized) + {developed_vp} (Developed) + {conquered_vp} (Conquered) + {research_vp} (Research)")
    print(f"Total VP = {total_vp}")

calculate_aog_vp()
<<<55>>>