def calculate_aog_vp():
    """
    Calculates and displays the total Victory Points (VP) for an Age of Galaxy game.
    """
    # --- Game State Variables ---
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
    # Terraforming, Advanced Terraforming, and Planetary Shields are all VP-scoring techs.
    researched_techs = 3

    # --- VP Calculation ---

    # Base VP
    vp_from_resources = (credits + productivity + discovery + influence) // 3
    vp_from_prestige = prestige
    vp_from_relics = relics
    vp_from_base_cruisers = cruisers
    vp_from_planets = (colonized_planets * 2) + (developed_planets * 1) + (conquered_planets * 1)
    vp_from_techs = researched_techs

    # Bonus VP from Alliances and Ideology
    vp_from_chaeilki = researched_techs // 2  # +1 VP per 2 techs
    vp_from_humans = developed_planets // 2    # +1 VP per 2 developed planets
    vp_from_usud = colonized_planets           # +1 VP per colonized planet
    vp_from_legarchaea = cruisers // 2         # +1 VP per 2 cruisers

    # Total VP
    total_vp = (
        vp_from_resources +
        vp_from_prestige +
        vp_from_relics +
        vp_from_base_cruisers +
        vp_from_planets +
        vp_from_techs +
        vp_from_chaeilki +
        vp_from_humans +
        vp_from_usud +
        vp_from_legarchaea
    )

    # --- Output ---
    print("Calculating Total Victory Points (VP)...")
    print("-" * 40)
    print("Final VP Equation:")
    print(
        f"{vp_from_resources} (from Resources) + "
        f"{vp_from_prestige} (from Prestige) + "
        f"{vp_from_relics} (from Relics) + "
        f"{vp_from_base_cruisers} (from Cruisers) + "
        f"{vp_from_planets} (from Planets) + "
        f"{vp_from_techs} (from Techs) + "
        f"{vp_from_chaeilki} (from Chaeilki Alliance) + "
        f"{vp_from_humans} (from Human Alliance) + "
        f"{vp_from_usud} (from Us'ud Alliance) + "
        f"{vp_from_legarchaea} (from Legarchaea Ideology)"
    )
    print("-" * 40)
    print(f"= {total_vp} Total VP")


if __name__ == "__main__":
    calculate_aog_vp()
<<<35>>>