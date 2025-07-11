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
total_vp = 0
vp_sources = []

# 1. Alliance VPs
# Chaeilki: 1 VP per planet
total_planets = colonized_planets + developed_planets + conquered_planets
chaeilki_vp = total_planets * 1
total_vp += chaeilki_vp
vp_sources.append(("Chaeilki Alliance (1 VP/planet)", chaeilki_vp))

# Humans: 2 VP per relic
humans_vp = relics * 2
total_vp += humans_vp
vp_sources.append(("Human Alliance (2 VP/relic)", humans_vp))

# Us'ud: 1 VP per research
usud_vp = researched_techs_count * 1
total_vp += usud_vp
vp_sources.append(("Us'ud Alliance (1 VP/research)", usud_vp))

# 2. Ideology VPs
# Legarchaea: 1 VP per influence
legarchaea_vp = influence * 1
total_vp += legarchaea_vp
vp_sources.append(("Legarchaea Ideology (1 VP/influence)", legarchaea_vp))

# 3. General VPs
# Colonized Planets: 1 VP each
colonized_vp = colonized_planets * 1
total_vp += colonized_vp
vp_sources.append(("Colonized Planets (1 VP each)", colonized_vp))

# Developed Planets: 2 VP each
developed_vp = developed_planets * 2
total_vp += developed_vp
vp_sources.append(("Developed Planets (2 VP each)", developed_vp))

# Conquered Planets: 3 VP each
conquered_vp = conquered_planets * 3
total_vp += conquered_vp
vp_sources.append(("Conquered Planets (3 VP each)", conquered_vp))

# Prestige: 1 VP per 5 prestige
prestige_vp = prestige // 5
total_vp += prestige_vp
vp_sources.append(("Prestige (1 VP/5 prestige)", prestige_vp))

# Relics: 1 VP each
relics_vp = relics * 1
total_vp += relics_vp
vp_sources.append(("Relics (1 VP each)", relics_vp))

# --- Final Output ---
print("Calculating Total Victory Points (VP):")
for reason, points in vp_sources:
    print(f"- {reason}: {points} VP")

print("\nFinal VP Equation:")
equation = " + ".join([str(points) for reason, points in vp_sources])
print(f"{equation} = {total_vp}")

print("\nTotal VP:", total_vp)
<<<26>>>