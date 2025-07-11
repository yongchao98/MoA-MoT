def calculate_age_of_galaxy_vp():
    """
    Calculates and displays the total Victory Points (VP) for a game of Age of Galaxy.
    """
    # --- Player's Assets ---
    prestige = 5
    relics = 2
    colonized_planets = 3
    developed_planets = 1
    conquered_planets = 1
    # Player has "Terraforming" and "Advanced Terraforming"
    terraforming_techs = 2

    # --- VP Calculation ---

    # 1. VP from Prestige
    vp_from_prestige = prestige

    # 2. VP from Base Planet Scores
    vp_from_colonized = colonized_planets * 1
    vp_from_developed_base = developed_planets * 3
    vp_from_conquered = conquered_planets * 2

    # 3. VP from Ideology (Legarchaea)
    # Rule: 2 VP per relic
    vp_from_legarchaea = relics * 2

    # 4. VP from Alliance Powers
    # Chaeilki: 1 VP per planet
    total_planets = colonized_planets + developed_planets + conquered_planets
    vp_from_chaeilki = total_planets * 1
    
    # Humans: 1 VP per developed planet
    vp_from_humans = developed_planets * 1
    
    # Us'ud: 1 VP per terraforming tech
    vp_from_usud = terraforming_techs * 1

    # --- Total VP ---
    total_vp = (
        vp_from_prestige + 
        vp_from_colonized + 
        vp_from_developed_base + 
        vp_from_conquered +
        vp_from_legarchaea +
        vp_from_chaeilki +
        vp_from_humans +
        vp_from_usud
    )

    # --- Print Detailed Breakdown ---
    print("Calculating Victory Points...")
    print("-" * 25)
    print(f"Prestige: {prestige} Prestige points = {vp_from_prestige} VP")
    print(f"Base Planet Score: {vp_from_colonized+vp_from_developed_base+vp_from_conquered} VP")
    print(f"  - {colonized_planets} Colonized Planet(s) * 1 VP = {vp_from_colonized} VP")
    print(f"  - {developed_planets} Developed Planet(s) * 3 VP = {vp_from_developed_base} VP")
    print(f"  - {conquered_planets} Conquered Planet(s) * 2 VP = {vp_from_conquered} VP")
    print(f"Ideology (Legarchaea): {relics} Relic(s) * 2 VP = {vp_from_legarchaea} VP")
    print("Alliance Score:")
    print(f"  - Chaeilki Bonus: {total_planets} total planets * 1 VP = {vp_from_chaeilki} VP")
    print(f"  - Human Bonus: {developed_planets} developed planet(s) * 1 VP = {vp_from_humans} VP")
    print(f"  - Us'ud Bonus: {terraforming_techs} terraforming techs * 1 VP = {vp_from_usud} VP")
    print("-" * 25)

    # --- Print Final Equation ---
    print("Final VP Equation:")
    final_equation = (
        f"{vp_from_prestige} (Prestige) + "
        f"{vp_from_colonized} (Colonized) + {vp_from_developed_base} (Developed) + {vp_from_conquered} (Conquered) + "
        f"{vp_from_legarchaea} (Legarchaea) + "
        f"{vp_from_chaeilki} (Chaeilki) + {vp_from_humans} (Humans) + {vp_from_usud} (Us'ud) = {total_vp}"
    )
    print(final_equation)
    print(f"\nTotal VP: {total_vp}")


if __name__ == "__main__":
    calculate_age_of_galaxy_vp()
<<<25>>>