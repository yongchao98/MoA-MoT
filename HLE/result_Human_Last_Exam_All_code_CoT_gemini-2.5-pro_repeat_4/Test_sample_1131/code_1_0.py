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
researched_terraforming = True
researched_advanced_terraforming = True
researched_planetary_shields = True

# --- VP Calculation ---

# 1. Alliance VPs
chaeilki_vp = relics
humans_vp = colonized_planets
usud_vp = developed_planets
alliance_vp = chaeilki_vp + humans_vp + usud_vp

# 2. Ideology VPs
legarchaea_vp = conquered_planets

# 3. General VPs
credits_vp = credits // 5
productivity_vp = productivity // 5
discovery_vp = discovery // 5
influence_vp = influence // 5
prestige_vp = prestige
relics_vp = relics
cruisers_vp = cruisers // 5

# 4. Planet VPs
colonized_planets_vp = colonized_planets * 1
developed_planets_vp = developed_planets * 2
conquered_planets_vp = conquered_planets * 1
planet_vp = colonized_planets_vp + developed_planets_vp + conquered_planets_vp

# 5. Research VPs
research_vp = 0
if researched_terraforming:
    research_vp += 1
if researched_advanced_terraforming:
    research_vp += 2
if researched_planetary_shields:
    research_vp += 1

# --- Total VP Calculation ---
total_vp = (
    alliance_vp +
    legarchaea_vp +
    credits_vp +
    productivity_vp +
    discovery_vp +
    influence_vp +
    prestige_vp +
    relics_vp +
    cruisers_vp +
    planet_vp +
    research_vp
)

# --- Output the breakdown ---
print("VP Calculation Breakdown:")
print(f"Alliance (Chaeilki from Relics): {chaeilki_vp}")
print(f"Alliance (Humans from Colonized Planets): {humans_vp}")
print(f"Alliance (Us'ud from Developed Planets): {usud_vp}")
print(f"Ideology (Legarchaea from Conquered Planets): {legarchaea_vp}")
print(f"Credits ({credits} // 5): {credits_vp}")
print(f"Productivity ({productivity} // 5): {productivity_vp}")
print(f"Discovery ({discovery} // 5): {discovery_vp}")
print(f"Influence ({influence} // 5): {influence_vp}")
print(f"Prestige: {prestige_vp}")
print(f"Relics: {relics_vp}")
print(f"Cruisers ({cruisers} // 5): {cruisers_vp}")
print(f"Colonized Planets ({colonized_planets} * 1): {colonized_planets_vp}")
print(f"Developed Planets ({developed_planets} * 2): {developed_planets_vp}")
print(f"Conquered Planets ({conquered_planets} * 1): {conquered_planets_vp}")
print(f"Research (Terraforming + Advanced Terraforming + Planetary Shields): {research_vp}")
print("-" * 25)
print("Final Equation:")
print(f"{chaeilki_vp} + {humans_vp} + {usud_vp} + {legarchaea_vp} + {credits_vp} + {productivity_vp} + {discovery_vp} + {influence_vp} + {prestige_vp} + {relics_vp} + {cruisers_vp} + {colonized_planets_vp} + {developed_planets_vp} + {conquered_planets_vp} + {research_vp} = {total_vp}")
print(f"\nTotal VP: {total_vp}")
