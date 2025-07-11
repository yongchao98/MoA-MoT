def calculate_aog_vp():
    """
    Calculates the total Victory Points (VP) for a game of Age of Galaxy
    based on the provided game state.
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

    # --- VP Calculation ---

    # 1. Alliance VP
    chaeilki_vp = relics * 3
    humans_vp = (colonized_planets * 2) + (developed_planets * 3)
    usud_vp = cruisers * 1

    # 2. Ideology VP (Legarchaea)
    legarchaea_planets = colonized_planets + developed_planets + conquered_planets
    legarchaea_vp = legarchaea_planets * 2

    # 3. Standard VP
    credits_vp = credits // 5
    productivity_vp = productivity // 5
    discovery_vp = discovery // 5
    influence_vp = influence // 5
    prestige_vp = prestige # 1 VP per prestige
    tech_vp = researched_techs # 1 VP per tech
    conquered_vp = conquered_planets * 2

    # Total VP
    total_vp = (
        chaeilki_vp + humans_vp + usud_vp + legarchaea_vp +
        credits_vp + productivity_vp + discovery_vp + influence_vp +
        prestige_vp + tech_vp + conquered_vp
    )

    # --- Outputting the result ---
    # We will build a string that shows the full equation
    # as requested, displaying each component.
    calculation_parts = [
        str(chaeilki_vp),      # Chaeilki Alliance (Relics)
        str(humans_vp),        # Human Alliance (Planets)
        str(usud_vp),          # Us'ud Alliance (Cruisers)
        str(legarchaea_vp),    # Legarchaea Ideology (Planets)
        str(prestige_vp),      # Prestige
        str(tech_vp),          # Researched Techs
        str(conquered_vp),     # Conquered Planets
        str(credits_vp),       # Credits
        str(productivity_vp),  # Productivity
        str(discovery_vp),     # Discovery
        str(influence_vp)      # Influence
    ]

    # Filter out zero-VP parts to keep the equation clean
    # For this specific case, discovery_vp and influence_vp are 0
    non_zero_parts = [part for part in calculation_parts if part != '0']

    equation = " + ".join(non_zero_parts)

    print(f"Total VP Breakdown:")
    print(f"Alliance (Chaeilki, Humans, Us'ud): {chaeilki_vp} + {humans_vp} + {usud_vp}")
    print(f"Ideology (Legarchaea): {legarchaea_vp}")
    print(f"Standard (Prestige, Tech, Conquests, Resources): {prestige_vp} + {tech_vp} + {conquered_vp} + {credits_vp} + {productivity_vp} + {discovery_vp} + {influence_vp}")
    print("-" * 20)
    print(f"Final Equation: {equation} = {total_vp}")
    print(f"Total VP: {total_vp}")

calculate_aog_vp()
<<<42>>>