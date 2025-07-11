def calculate_aog_vp():
    """
    Calculates and prints the total Victory Points (VP) for a game of Age of Galaxy.
    """
    # Player's assets
    credits = 7
    productivity = 5
    discovery = 4
    influence = 3
    prestige = 5
    relics = 2
    cruisers = 5
    colonized_planets = 3
    developed_planets = 1
    conquered_planets = 1
    researched_techs = 3 # Terraforming, Advanced Terraforming, Planetary Shields

    # --- VP Calculations ---

    # 1. Alliance VP
    # Chaeilki: 1 VP per 2 Cruisers
    chaeilki_vp = cruisers // 2
    # Humans: 1 VP per 2 Colonized Planets
    humans_vp = colonized_planets // 2
    # Us'ud: 1 VP per 2 Relics
    usud_vp = relics // 2
    total_alliance_vp = chaeilki_vp + humans_vp + usud_vp

    # 2. Ideology VP
    # Legarchaea: 1 VP per 2 Prestige
    legarchaea_vp = prestige // 2

    # 3. Resource VP
    # 1 VP per 7 credits, 5 productivity, 4 discovery, 3 influence
    credits_vp = credits // 7
    productivity_vp = productivity // 5
    discovery_vp = discovery // 4
    influence_vp = influence // 3
    total_resource_vp = credits_vp + productivity_vp + discovery_vp + influence_vp

    # 4. Planet VP
    # Colonized: 1 VP, Developed: 2 VP, Conquered: 3 VP
    colonized_vp = colonized_planets * 1
    developed_vp = developed_planets * 2
    conquered_vp = conquered_planets * 3
    total_planet_vp = colonized_vp + developed_vp + conquered_vp

    # 5. Research VP
    # 1 VP per researched tech
    research_vp = researched_techs * 1

    # --- Total VP ---
    total_vp = total_alliance_vp + legarchaea_vp + total_resource_vp + total_planet_vp + research_vp

    # --- Print the breakdown ---
    print("VP Calculation Breakdown:")
    print(f"Alliance VP: {chaeilki_vp} (Chaeilki) + {humans_vp} (Humans) + {usud_vp} (Us'ud) = {total_alliance_vp}")
    print(f"Ideology VP (Legarchaea): {legarchaea_vp}")
    print(f"Resource VP: {credits_vp} (Credits) + {productivity_vp} (Productivity) + {discovery_vp} (Discovery) + {influence_vp} (Influence) = {total_resource_vp}")
    print(f"Planet VP: {colonized_vp} (Colonized) + {developed_vp} (Developed) + {conquered_vp} (Conquered) = {total_planet_vp}")
    print(f"Research VP: {research_vp}")
    print("-" * 20)
    print(f"Total VP = {total_alliance_vp} + {legarchaea_vp} + {total_resource_vp} + {total_planet_vp} + {research_vp} = {total_vp}")

calculate_aog_vp()
<<<21>>>