# Player's game state
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

# 1. Alliance VPs
# Chaeilki: 1 VP per 2 Productivity
chaeilki_vp = productivity // 2
# Humans: 1 VP per 3 colonized planets
humans_vp = colonized_planets // 3
# Us'ud: 1 VP per 2 Prestige
usud_vp = prestige // 2

# 2. Ideology VPs
# Legarchaea (Major): 2 VP per Relic
legarchaea_vp = relics * 2

# 3. Standard VPs
# Note: Base game VP for resources is 1 VP per 10. These will be 0.
credits_vp = credits // 10
productivity_vp = productivity // 10
discovery_vp = discovery // 10
influence_vp = influence // 10

# Other standard VPs
prestige_vp = prestige
relics_vp = relics
cruisers_vp = cruisers // 5
colonized_planets_vp = colonized_planets * 2
developed_planets_vp = developed_planets * 3
conquered_planets_vp = conquered_planets * 2
researched_techs_vp = researched_techs * 1

# --- Total VP Calculation ---
total_vp = (
    chaeilki_vp + humans_vp + usud_vp + legarchaea_vp +
    credits_vp + productivity_vp + discovery_vp + influence_vp +
    prestige_vp + relics_vp + cruisers_vp +
    colonized_planets_vp + developed_planets_vp + conquered_planets_vp +
    researched_techs_vp
)

# --- Output the result ---
print("Calculating Victory Points...\n")

print("Alliance VPs:")
print(f"  Chaeilki (5 Productivity // 2): {chaeilki_vp} VP")
print(f"  Humans (3 Colonized Planets // 3): {humans_vp} VP")
print(f"  Us'ud (5 Prestige // 2): {usud_vp} VP")
print("\nIdeology VPs:")
print(f"  Legarchaea (2 Relics * 2): {legarchaea_vp} VP")
print("\nStandard VPs:")
print(f"  Resources (Credits, Productivity, etc.): {credits_vp + productivity_vp + discovery_vp + influence_vp} VP")
print(f"  Prestige: {prestige_vp} VP")
print(f"  Relics: {relics_vp} VP")
print(f"  Cruisers (5 Cruisers // 5): {cruisers_vp} VP")
print(f"  Colonized Planets (3 * 2): {colonized_planets_vp} VP")
print(f"  Developed Planets (1 * 3): {developed_planets_vp} VP")
print(f"  Conquered Planets (1 * 2): {conquered_planets_vp} VP")
print(f"  Researched Techs (3 * 1): {researched_techs_vp} VP")
print("\n---")
print("Final VP Calculation:")
print(f"{chaeilki_vp} (Chaeilki) + {humans_vp} (Humans) + {usud_vp} (Us'ud) + {legarchaea_vp} (Legarchaea) + {prestige_vp} (Prestige) + {relics_vp} (Relics) + {cruisers_vp} (Cruisers) + {colonized_planets_vp} (Colonized) + {developed_planets_vp} (Developed) + {conquered_planets_vp} (Conquered) + {researched_techs_vp} (Techs) = {total_vp}")
print(f"\nTotal VP: {total_vp}")
