def calculate_age_of_galaxy_vp():
    """
    Calculates the total Victory Points (VP) for a player in Age of Galaxy
    based on their end-game state.
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

    total_vp = 0
    calculation_steps = []

    # --- Base VP Calculation ---
    # Influence VP
    vp_influence = influence * 1
    total_vp += vp_influence
    calculation_steps.append(f"{vp_influence} (from {influence} Influence)")

    # Prestige VP
    vp_prestige_base = prestige * 1
    total_vp += vp_prestige_base
    calculation_steps.append(f"{vp_prestige_base} (from {prestige} Prestige)")
    
    # Relics VP
    vp_relics_base = relics * 3
    total_vp += vp_relics_base
    calculation_steps.append(f"{vp_relics_base} (from {relics} Relics * 3)")
    
    # Discovery VP
    vp_discovery = discovery // 2
    total_vp += vp_discovery
    calculation_steps.append(f"{vp_discovery} (from {discovery} Discovery // 2)")
    
    # Planet VP
    vp_colonized_base = colonized_planets * 1
    total_vp += vp_colonized_base
    calculation_steps.append(f"{vp_colonized_base} (from {colonized_planets} Colonized Planets * 1)")
    
    vp_developed_base = developed_planets * 2
    total_vp += vp_developed_base
    calculation_steps.append(f"{vp_developed_base} (from {developed_planets} Developed Planet * 2)")

    vp_conquered_base = conquered_planets * 1
    total_vp += vp_conquered_base
    calculation_steps.append(f"{vp_conquered_base} (from {conquered_planets} Conquered Planet * 1)")

    # --- Alliance VP Calculation ---
    # Chaeilki Alliance
    vp_chaeilki_relics = relics * 1
    total_vp += vp_chaeilki_relics
    calculation_steps.append(f"{vp_chaeilki_relics} (from {relics} Relics [Chaeilki])")

    vp_chaeilki_prestige = prestige * 1
    total_vp += vp_chaeilki_prestige
    calculation_steps.append(f"{vp_chaeilki_prestige} (from {prestige} Prestige [Chaeilki])")

    # Humans Alliance
    vp_humans_colonized = colonized_planets * 1
    total_vp += vp_humans_colonized
    calculation_steps.append(f"{vp_humans_colonized} (from {colonized_planets} Colonized Planets [Humans])")

    vp_humans_developed = developed_planets * 1
    total_vp += vp_humans_developed
    calculation_steps.append(f"{vp_humans_developed} (from {developed_planets} Developed Planet [Humans])")
    
    # Us'ud Alliance
    vp_usud_credits = credits // 3
    total_vp += vp_usud_credits
    calculation_steps.append(f"{vp_usud_credits} (from {credits} Credits // 3 [Us'ud])")
    
    vp_usud_productivity = productivity // 3
    total_vp += vp_usud_productivity
    calculation_steps.append(f"{vp_usud_productivity} (from {productivity} Productivity // 3 [Us'ud])")

    # --- Ideology VP Calculation ---
    # Legarchaea Ideology
    vp_legarchaea_conquered = conquered_planets * 1
    total_vp += vp_legarchaea_conquered
    calculation_steps.append(f"{vp_legarchaea_conquered} (from {conquered_planets} Conquered Planet [Legarchaea])")
    
    vp_legarchaea_cruisers = cruisers * 1
    total_vp += vp_legarchaea_cruisers
    calculation_steps.append(f"{vp_legarchaea_cruisers} (from {cruisers} Cruisers [Legarchaea])")
    
    # --- Final Output ---
    print("VP Calculation Breakdown:")
    final_equation = " + ".join([step.split(" ")[0] for step in calculation_steps])
    print(f"{final_equation} = {total_vp}")
    
    print("\nDetails:")
    for step in calculation_steps:
        print(f"- {step}")
        
    print(f"\nTotal VP: {total_vp}")

calculate_age_of_galaxy_vp()
<<<42>>>