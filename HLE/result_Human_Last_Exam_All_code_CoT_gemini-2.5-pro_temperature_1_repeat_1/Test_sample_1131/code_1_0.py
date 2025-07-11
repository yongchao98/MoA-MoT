# Game state variables
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
# Researched technologies do not grant VPs themselves, but Planetary Shields has an effect.

# --- VP Calculation ---

# 1. Alliance VPs
chaeilki_vp = colonized_planets * 2  # 2 VP per colonized planet
humans_vp = developed_planets * 3    # 3 VP per developed planet
usud_vp = relics * 2                 # 2 VP per relic
alliance_total_vp = chaeilki_vp + humans_vp + usud_vp

# 2. Ideology VPs (Legarchaea)
legarchaea_vp = cruisers * 2 # 2 VP per cruiser

# 3. Planet VPs
colonized_planet_vp = colonized_planets * 1
developed_planet_vp = developed_planets * 2
conquered_planet_vp = conquered_planets * 1
planet_total_vp = colonized_planet_vp + developed_planet_vp + conquered_planet_vp

# 4. Resource VPs
standard_resources_vp = (credits + productivity + discovery + influence) // 5
prestige_vp = prestige * 1
resource_total_vp = standard_resources_vp + prestige_vp

# 5. Technology VPs (Planetary Shields)
planetary_shields_vp = developed_planets * 1 # 1 VP per developed planet

# --- Total VP Calculation ---
total_vp = alliance_total_vp + legarchaea_vp + planet_total_vp + resource_total_vp + planetary_shields_vp

# --- Output ---
print("Calculating Total Victory Points (VP):")
print("\n--- Alliance VP ---")
print(f"Chaeilki: {chaeilki_vp} VP ({colonized_planets} colonized planets * 2 VP)")
print(f"Humans: {humans_vp} VP ({developed_planets} developed planet * 3 VP)")
print(f"Us'ud: {usud_vp} VP ({relics} relics * 2 VP)")
print(f"Subtotal: {alliance_total_vp} VP")

print("\n--- Ideology VP (Legarchaea) ---")
print(f"Legarchaea: {legarchaea_vp} VP ({cruisers} cruisers * 2 VP)")
print(f"Subtotal: {legarchaea_vp} VP")

print("\n--- Planet VP ---")
print(f"Colonized Planets: {colonized_planet_vp} VP ({colonized_planets} * 1 VP)")
print(f"Developed Planets: {developed_planet_vp} VP ({developed_planets} * 2 VP)")
print(f"Conquered Planets: {conquered_planet_vp} VP ({conquered_planets} * 1 VP)")
print(f"Subtotal: {planet_total_vp} VP")

print("\n--- Resource VP ---")
print(f"Standard Resources: {standard_resources_vp} VP (({credits} + {productivity} + {discovery} + {influence}) / 5)")
print(f"Prestige: {prestige_vp} VP ({prestige} * 1 VP)")
print(f"Subtotal: {resource_total_vp} VP")

print("\n--- Technology VP (Planetary Shields) ---")
print(f"Planetary Shields: {planetary_shields_vp} VP ({developed_planets} developed planet * 1 VP)")
print(f"Subtotal: {planetary_shields_vp} VP")

print("\n--- Final Score ---")
print(f"Total VP = {alliance_total_vp} (Alliances) + {legarchaea_vp} (Ideology) + {planet_total_vp} (Planets) + {resource_total_vp} (Resources) + {planetary_shields_vp} (Tech)")
print(f"Total VP = {total_vp}")
