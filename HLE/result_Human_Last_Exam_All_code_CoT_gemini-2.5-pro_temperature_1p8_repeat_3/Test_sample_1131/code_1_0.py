def calculate_age_of_galaxy_vp():
    """
    Calculates the total Victory Points (VP) for a player in Age of Galaxy
    based on their alliance, ideology, and assets.
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

    # --- VP Calculation ---
    total_vp = 0
    calculation_parts = []
    
    print("Calculating Total Victory Points (VP)...")
    print("-" * 30)

    # 1. VP from Prestige (Us'ud bonus: +1 VP per prestige)
    prestige_vp = prestige * 1
    total_vp += prestige_vp
    calculation_parts.append(str(prestige_vp))
    print(f"VP from Prestige (Us'ud): {prestige} prestige * 1 VP/prestige = {prestige_vp} VP")

    # 2. VP from Relics (Chaeilki: +1 VP, Legarchaea: +2 VP)
    # Total bonus per relic is 1 + 2 = 3 VP
    relics_vp = relics * (1 + 2)
    total_vp += relics_vp
    calculation_parts.append(str(relics_vp))
    print(f"VP from Relics (Chaeilki + Legarchaea): {relics} relics * 3 VP/relic = {relics_vp} VP")

    # 3. VP from Productivity (Chaeilki bonus: +1 VP per productivity)
    productivity_vp = productivity * 1
    total_vp += productivity_vp
    calculation_parts.append(str(productivity_vp))
    print(f"VP from Productivity (Chaeilki): {productivity} productivity * 1 VP/prod = {productivity_vp} VP")

    # 4. VP from Discovery (Human bonus: +1 VP per discovery)
    discovery_vp = discovery * 1
    total_vp += discovery_vp
    calculation_parts.append(str(discovery_vp))
    print(f"VP from Discovery (Humans): {discovery} discovery * 1 VP/discovery = {discovery_vp} VP")

    # 5. VP from Planets
    # Human bonus: +1 VP per colonized planet.
    # It is assumed that "3 colonized planets" is the total, and the developed planet is one of these.
    colonized_vp = colonized_planets * 1
    total_vp += colonized_vp
    calculation_parts.append(str(colonized_vp))
    print(f"VP from Colonized Planets (Humans): {colonized_planets} planets * 1 VP/planet = {colonized_vp} VP")
    
    # Us'ud bonus: +1 VP per developed planet. This is an additional point.
    developed_vp = developed_planets * 1
    total_vp += developed_vp
    calculation_parts.append(str(developed_vp))
    print(f"VP from Developed Planet (Us'ud): {developed_planets} planet * 1 VP/planet = {developed_vp} VP")

    # Legarchaea penalty: -1 VP per conquered planet
    conquered_vp = conquered_planets * -1
    total_vp += conquered_vp
    calculation_parts.append(str(conquered_vp))
    print(f"VP from Conquered Planet (Legarchaea): {conquered_planets} planet * -1 VP/planet = {conquered_vp} VP")

    # 6. VP from other sources
    # Cruisers, researched techs, credits, and influence have no scoring rules provided.
    print("VP from Cruisers, Research, and other resources: 0 VP (no scoring rules provided)")

    print("-" * 30)
    final_equation = " + ".join(calculation_parts).replace("+ -", "- ")
    print(f"Final Calculation: {final_equation} = {total_vp}")
    print(f"\nThe player's total Victory Points are: {total_vp}")

# Run the calculation and print the result
calculate_age_of_galaxy_vp()