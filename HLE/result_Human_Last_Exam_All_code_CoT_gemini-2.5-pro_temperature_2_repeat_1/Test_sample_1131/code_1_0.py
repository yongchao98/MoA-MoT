import math

# Player's assets and game state
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
num_factions_in_alliance = 3
# Technologies do not provide direct VP in this calculation
techs_researched = ["Terraforming", "Advanced Terraforming", "Planetary Shields"]

# Initialize total VP and a list to hold the numbers for the final equation
total_vp = 0
vp_components = []

print("Calculating Total Victory Points (VP)...\n")

# --- Alliance & Faction VP ---
# Chaeilki: 1 VP per colonized planet
chaeilki_vp = 1 * colonized_planets
total_vp += chaeilki_vp
vp_components.append(chaeilki_vp)
print(f"Chaeilki Faction VP (1 per colonized planet): {chaeilki_vp}")

# Humans: 1 VP per faction in the alliance
humans_vp = 1 * num_factions_in_alliance
total_vp += humans_vp
vp_components.append(humans_vp)
print(f"Human Faction VP (1 per alliance member): {humans_vp}")

# Us'ud: 2 VP per relic
usud_vp = 2 * relics
total_vp += usud_vp
vp_components.append(usud_vp)
print(f"Us'ud Faction VP (2 per relic): {usud_vp}")

# Alliance of 3 Bonus
alliance_size_vp = 3
total_vp += alliance_size_vp
vp_components.append(alliance_size_vp)
print(f"Alliance Size Bonus (3 members): {alliance_size_vp}")

# --- Ideology VP ---
# Legarchaea: 2 VP per colonized planet
legarchaea_vp = 2 * colonized_planets
total_vp += legarchaea_vp
vp_components.append(legarchaea_vp)
print(f"Legarchaea Ideology VP (2 per colonized planet): {legarchaea_vp}")

# --- General VP ---
# Resources (CPDI): 1 VP per 10 total resources
total_resources = credits + productivity + discovery + influence
resource_vp = total_resources // 10
total_vp += resource_vp
vp_components.append(resource_vp)
print(f"Resource VP (1 per 10 CPDI): {resource_vp}")

# Prestige: 1 VP per prestige point
prestige_vp = prestige
total_vp += prestige_vp
vp_components.append(prestige_vp)
print(f"Prestige VP (1 per prestige): {prestige_vp}")

# Relics (Base Value): 1 VP per relic
relic_vp_base = relics
total_vp += relic_vp_base
vp_components.append(relic_vp_base)
print(f"Base Relic VP (1 per relic): {relic_vp_base}")

# Cruisers: 1 VP per 3 cruisers
cruiser_vp = cruisers // 3
total_vp += cruiser_vp
vp_components.append(cruiser_vp)
print(f"Cruiser VP (1 per 3 cruisers): {cruiser_vp}")

# Planets (Base Value)
colonized_planet_vp = 1 * colonized_planets
total_vp += colonized_planet_vp
vp_components.append(colonized_planet_vp)
print(f"Colonized Planet VP (1 per planet): {colonized_planet_vp}")

developed_planet_vp = 3 * developed_planets
total_vp += developed_planet_vp
vp_components.append(developed_planet_vp)
print(f"Developed Planet VP (3 per planet): {developed_planet_vp}")

conquered_planet_vp = 1 * conquered_planets
total_vp += conquered_planet_vp
vp_components.append(conquered_planet_vp)
print(f"Conquered Planet VP (1 per planet): {conquered_planet_vp}")


# --- Final Tally ---
print("\n" + "="*25)
# Using map to convert all numbers in the list to strings for the join function
final_equation_str = " + ".join(map(str, vp_components))
print(f"Final Equation: {final_equation_str}")
print(f"GRAND TOTAL: {total_vp} VP")
print("="*25)