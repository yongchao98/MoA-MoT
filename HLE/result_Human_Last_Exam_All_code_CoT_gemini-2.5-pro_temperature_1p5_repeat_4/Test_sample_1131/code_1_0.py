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
researched_techs_count = 3 # Terraforming, Advanced Terraforming, Planetary Shields

# --- VP Calculation ---

# 1. Alliance VP
# Chaeilki: 1 VP per Relic
chaeilki_vp = 1 * relics
# Humans: 1 VP per planet (colonized, developed, conquered)
total_planets = colonized_planets + developed_planets + conquered_planets
humans_vp = 1 * total_planets
# Us'ud: 2 VP per Cruiser
usud_vp = 2 * cruisers

# 2. Ideology VP
# Legarchaea: 1 VP per Prestige
legarchaea_vp = 1 * prestige

# 3. General VP
credits_vp = credits // 5
productivity_vp = productivity // 5
discovery_vp = discovery // 5
influence_vp = influence // 5
tech_vp = 1 * researched_techs_count

# --- Total VP ---
total_vp = (chaeilki_vp + humans_vp + usud_vp + 
            legarchaea_vp + 
            credits_vp + productivity_vp + discovery_vp + influence_vp + tech_vp)

# --- Output ---
print("Calculating Total Victory Points (VP):")
print("-" * 40)
print("Alliance VP:")
print(f"  Chaeilki (2 Relics * 1 VP/Relic) = {chaeilki_vp} VP")
print(f"  Humans ({total_planets} Planets * 1 VP/Planet) = {humans_vp} VP")
print(f"  Us'ud (5 Cruisers * 2 VP/Cruiser) = {usud_vp} VP")
print("\nIdeology VP:")
print(f"  Legarchaea (5 Prestige * 1 VP/Prestige) = {legarchaea_vp} VP")
print("\nGeneral VP:")
print(f"  Resources (7Cr, 5P, 4D, 3I) = {credits_vp} + {productivity_vp} + {discovery_vp} + {influence_vp} VP")
print(f"  Researched Techs (3 Techs * 1 VP/Tech) = {tech_vp} VP")
print("-" * 40)
print(f"Final Equation: {chaeilki_vp} + {humans_vp} + {usud_vp} + {legarchaea_vp} + {credits_vp} + {productivity_vp} + {discovery_vp} + {influence_vp} + {tech_vp} = {total_vp}")
print(f"\nTotal VP = {total_vp}")
