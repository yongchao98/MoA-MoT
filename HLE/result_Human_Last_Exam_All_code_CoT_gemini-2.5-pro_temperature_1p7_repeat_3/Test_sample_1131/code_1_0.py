def calculate_age_of_galaxy_vp():
    """
    Calculates the total Victory Points (VP) for a player in Age of Galaxy
    based on their end-game state.
    """
    # Player's Assets
    num_factions = 3  # Chaeilki, Humans, Us'ud
    prestige = 5
    relics = 2
    cruisers = 5
    colonized_planets = 3
    developed_planets = 1
    conquered_planets = 1

    # Total planets under control
    total_planets = colonized_planets + developed_planets + conquered_planets

    # VP Calculations based on game rules
    # Faction VPs
    chaeilki_vp = developed_planets * 2
    humans_vp = total_planets * 1
    usud_vp = (cruisers // 2) * 1

    # Ideology VP
    legarchaea_vp = num_factions * 2

    # Standard VPs
    prestige_vp = prestige * 1
    relics_vp = relics * 3

    # Technology VP
    planetary_shields_vp = total_planets * 1

    # Store VP sources for clear output
    vp_sources = {
        "Chaeilki (Developed Planets)": chaeilki_vp,
        "Humans (Total Planets)": humans_vp,
        "Us'ud (Cruisers)": usud_vp,
        "Legarchaea Ideology (Factions)": legarchaea_vp,
        "Prestige": prestige_vp,
        "Relics": relics_vp,
        "Planetary Shields Tech (Planets)": planetary_shields_vp,
    }

    # Calculate total VP
    total_vp = sum(vp_sources.values())

    # Build the final equation string
    equation_parts = []
    for source, value in vp_sources.items():
        equation_parts.append(str(value))

    equation_str = " + ".join(equation_parts)

    # Print the results
    print("Calculating Total Victory Points (VP):")
    for source, value in vp_sources.items():
        print(f"- {source}: {value} VP")

    print("\nFinal VP Calculation:")
    print(f"{equation_str} = {total_vp}")
    
    print(f"\nTotal VP = {total_vp}")
    print(f"<<<{total_vp}>>>")


if __name__ == "__main__":
    calculate_age_of_galaxy_vp()