def calculate_aog_vp():
    """
    Calculates and prints the total Victory Points (VP) for a player in Age of Galaxy.
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

    # This dictionary will store the description and calculated VP for each source.
    vp_sources = {}

    # --- Alliance & Ideology VP ---
    vp_sources["Chaeilki Alliance (1 per developed planet)"] = 1 * developed_planets
    vp_sources["Human Alliance (1 per colonized planet)"] = 1 * colonized_planets
    vp_sources["Us'ud Alliance (2 per relic)"] = 2 * relics
    vp_sources["Legarchaea Ideology (1 per cruiser)"] = 1 * cruisers
    
    # --- Standard VP from assets ---
    vp_sources["Prestige"] = prestige
    vp_sources["Colonized Planets (1 VP each)"] = 1 * colonized_planets
    vp_sources["Developed Planets (2 VP each)"] = 2 * developed_planets
    vp_sources["Conquered Planets (1 VP each)"] = 1 * conquered_planets
    vp_sources["Relics (1 VP each)"] = 1 * relics

    # --- VP from resources ---
    vp_sources["Credits (1 per 7)"] = credits // 7
    total_resources = productivity + discovery + influence
    vp_sources["Prod/Disc/Inf (1 per 7 combined)"] = total_resources // 7

    # --- VP from researched techs ---
    vp_sources["Terraforming Tech"] = 1
    vp_sources["Advanced Terraforming Tech"] = 2
    vp_sources["Planetary Shields Tech"] = 1
    
    # --- Calculate Total and build the equation string ---
    total_vp = 0
    equation_parts = []
    
    print("Calculating VP from all sources:\n")
    for source, value in vp_sources.items():
        print(f"- {source}: {value} VP")
        total_vp += value
        equation_parts.append(str(value))
        
    final_equation = " + ".join(equation_parts)

    print("\n--- Total VP Calculation ---")
    print(f"Final Equation: {final_equation}")
    print(f"Total VP = {total_vp}")

calculate_aog_vp()
<<<32>>>