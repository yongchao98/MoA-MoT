import math

# --- Player's Assets ---
prestige = 5
relics = 2
discovery = 4
colonized_planets = 3
developed_planets = 1
conquered_planets = 1
alliance = ["Chaeilki", "Humans", "Us'ud"]

# --- VP Calculation ---
# A list to store the individual VP amounts for the final equation
vp_components = []
total_vp = 0

print("Calculating Victory Points...\n")

# 1. Ideology: Legarchaea
# 2 VP per faction in the alliance
ideology_faction_vp = len(alliance) * 2
vp_components.append(ideology_faction_vp)
total_vp += ideology_faction_vp
print(f"Legarchaea (Factions): {len(alliance)} factions * 2 VP = {ideology_faction_vp} VP")

# 3 VP per developed planet
ideology_developed_vp = developed_planets * 3
vp_components.append(ideology_developed_vp)
total_vp += ideology_developed_vp
print(f"Legarchaea (Developed Planets): {developed_planets} developed * 3 VP = {ideology_developed_vp} VP")

# 2. Standard Planet VP
# 1 VP per colonized planet
planet_colonized_vp = colonized_planets * 1
vp_components.append(planet_colonized_vp)
total_vp += planet_colonized_vp
print(f"Standard (Colonized Planets): {colonized_planets} colonized * 1 VP = {planet_colonized_vp} VP")

# 2 VP per developed planet
planet_developed_vp = developed_planets * 2
vp_components.append(planet_developed_vp)
total_vp += planet_developed_vp
print(f"Standard (Developed Planets): {developed_planets} developed * 2 VP = {planet_developed_vp} VP")

# 1 VP per conquered planet
planet_conquered_vp = conquered_planets * 1
vp_components.append(planet_conquered_vp)
total_vp += planet_conquered_vp
print(f"Standard (Conquered Planets): {conquered_planets} conquered * 1 VP = {planet_conquered_vp} VP")

# 3. Prestige VP
# 1 VP per prestige
prestige_vp = prestige * 1
vp_components.append(prestige_vp)
total_vp += prestige_vp
print(f"Prestige: {prestige} prestige * 1 VP = {prestige_vp} VP")

# 4. Relic VP
# 3 VP per relic
relic_vp = relics * 3
vp_components.append(relic_vp)
total_vp += relic_vp
print(f"Relics: {relics} relics * 3 VP = {relic_vp} VP")

# 5. Alliance Bonuses
print("\n--- Alliance Bonuses ---")
# Chaeilki: 1 VP per 2 discovery
chaeilki_vp = math.floor(discovery / 2) * 1
vp_components.append(chaeilki_vp)
total_vp += chaeilki_vp
print(f"Chaeilki Bonus (Discovery): floor({discovery} / 2) * 1 VP = {chaeilki_vp} VP")

# Humans: 1 VP per colonized planet
humans_vp = colonized_planets * 1
vp_components.append(humans_vp)
total_vp += humans_vp
print(f"Human Bonus (Colonized): {colonized_planets} colonized * 1 VP = {humans_vp} VP")

# Us'ud: 1 VP per 2 relics
usud_vp = math.floor(relics / 2) * 1
vp_components.append(usud_vp)
total_vp += usud_vp
print(f"Us'ud Bonus (Relics): floor({relics} / 2) * 1 VP = {usud_vp} VP")

# --- Final Score ---
# Create the equation string
equation_string = " + ".join(map(str, vp_components))

print("\n--- Total VP ---")
print(f"Final Calculation: {equation_string} = {total_vp}")
print(f"The player's total VP is: {total_vp}")
