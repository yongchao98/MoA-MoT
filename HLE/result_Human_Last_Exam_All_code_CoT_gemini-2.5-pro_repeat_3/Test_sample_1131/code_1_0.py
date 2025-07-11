def calculate_aog_vp():
    """
    Calculates and prints the total Victory Points (VP) for an Age of Galaxy game session.
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
    num_techs = 3 # Terraforming, Advanced Terraforming, Planetary Shields
    num_alliances = 3 # Chaeilki, Humans, Us'ud

    # VP Calculations

    # 1. Prestige VP: 1 VP per prestige point
    prestige_vp = prestige * 1

    # 2. Relic VP: 3 VP per relic
    relic_vp = relics * 3

    # 3. Planet VP: 1 VP per colonized/conquered, 3 VP per developed
    planet_vp = (colonized_planets * 1) + (developed_planets * 3) + (conquered_planets * 1)

    # 4. Resource VP: 1 VP for each set of 3 resources
    total_resources = credits + productivity + discovery + influence
    resource_vp = total_resources // 3

    # 5. Technology VP: 1 VP per researched technology
    tech_vp = num_techs * 1

    # 6. Legarchaea (Ideology) VP: 3 VP for each alliance
    legarchaea_vp = num_alliances * 3

    # 7. Chaeilki (Alliance) VP: 2 VP for each developed planet
    chaeilki_vp = developed_planets * 2

    # 8. Humans (Alliance) VP: 1 VP for each colonized and conquered planet
    humans_vp = (colonized_planets + conquered_planets) * 1

    # 9. Us'ud (Alliance) VP: 1 VP for every 2 cruisers
    usud_vp = cruisers // 2

    # Total VP
    total_vp = (prestige_vp + relic_vp + planet_vp + resource_vp + tech_vp +
                legarchaea_vp + chaeilki_vp + humans_vp + usud_vp)

    # Print the detailed equation
    equation = [
        f"{prestige_vp} (Prestige)",
        f"{relic_vp} (Relics)",
        f"{planet_vp} (Planets)",
        f"{resource_vp} (Resources)",
        f"{tech_vp} (Technologies)",
        f"{legarchaea_vp} (Legarchaea Ideology)",
        f"{chaeilki_vp} (Chaeilki Alliance)",
        f"{humans_vp} (Humans Alliance)",
        f"{usud_vp} (Us'ud Alliance)"
    ]

    print("VP Calculation Breakdown:")
    print(f"Total VP = {' + '.join(equation)} = {total_vp}")
    
# Execute the calculation
calculate_aog_vp()