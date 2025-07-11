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

# 1. Alliance VPs
# 3-member alliance = 9 VP
alliance_size_vp = 9
# Chaeilki: +1 VP per planet (3 colonized + 1 developed + 1 conquered = 5 planets)
chaeilki_vp = colonized_planets + developed_planets + conquered_planets
# Humans: +1 VP per relic
humans_vp = relics
# Us'ud: +1 VP per prestige
usud_vp = prestige
total_alliance_vp = alliance_size_vp + chaeilki_vp + humans_vp + usud_vp

# 2. Ideology VPs (Legarchaea)
# +2 VP per developed planet
ideology_vp = 2 * developed_planets

# 3. Planet VPs
# Colonized: 1 VP each, Developed: 3 VP each, Conquered: 2 VP each
planet_vp = (colonized_planets * 1) + (developed_planets * 3) + (conquered_planets * 2)

# 4. Resource VPs (lowest of Credits, Productivity, Discovery, Influence)
resource_vp = min(credits, productivity, discovery, influence)

# 5. Prestige & Relic VPs
prestige_vp = prestige
relic_vp = relics

# 6. Military VPs
# Cruisers: 1 VP each
military_vp = cruisers * 1

# 7. Research VPs
# Terraforming (1) + Advanced Terraforming (2) + Planetary Shields (2)
research_vp = 1 + 2 + 2

# 8. Total VPs
total_vp = total_alliance_vp + ideology_vp + planet_vp + resource_vp + prestige_vp + relic_vp + military_vp + research_vp

# --- Output the result ---
print("Calculating Total Victory Points (VP):")
print(f"Alliance VP: {alliance_size_vp} (3 members) + {chaeilki_vp} (Chaeilki) + {humans_vp} (Humans) + {usud_vp} (Us'ud) = {total_alliance_vp}")
print(f"Ideology VP (Legarchaea): {ideology_vp}")
print(f"Planet VP: {colonized_planets * 1} (Colonized) + {developed_planets * 3} (Developed) + {conquered_planets * 2} (Conquered) = {planet_vp}")
print(f"Resource VP: {resource_vp}")
print(f"Prestige VP: {prestige_vp}")
print(f"Relic VP: {relic_vp}")
print(f"Military VP: {military_vp}")
print(f"Research VP: {research_vp}")
print("-" * 20)
print(f"Final Equation: {total_alliance_vp} + {ideology_vp} + {planet_vp} + {resource_vp} + {prestige_vp} + {relic_vp} + {military_vp} + {research_vp} = {total_vp}")
print(f"Total VP: {total_vp}")
