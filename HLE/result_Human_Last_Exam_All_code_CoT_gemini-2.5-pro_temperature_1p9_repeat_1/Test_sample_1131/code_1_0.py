import math

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
techs = ["Terraforming", "Advanced Terraforming", "Planetary Shields"]

# VP Calculation based on standard rules

# 1. Base VP from Planets
# 1 VP per Colonized, 3 per Developed, 2 per Conquered
vp_colonized = colonized_planets * 1
vp_developed = developed_planets * 3
vp_conquered = conquered_planets * 2
vp_from_planets_base = vp_colonized + vp_developed + vp_conquered

# 2. VP from Technology
# 2 VP for Advanced Terraforming, 2 for Planetary Shields
vp_from_tech = 0
if "Advanced Terraforming" in techs:
    vp_from_tech += 2
if "Planetary Shields" in techs:
    vp_from_tech += 2

# 3. VP from Military
# 1 VP per 2 cruisers
vp_from_military = math.floor(cruisers / 2)

# 4. Base VP from Relics
# 3 VP per Relic
vp_from_relics_base = relics * 3

# 5. VP from Prestige
# 1 VP per 5 Prestige
vp_from_prestige = math.floor(prestige / 5)

# 6. Alliance and Ideology Bonus VP
# Chaeilki: +1 VP per planet
total_planets = colonized_planets + developed_planets + conquered_planets
vp_bonus_chaeilki = total_planets * 1

# Humans: +1 VP per resource type (credits, productivity, discovery, influence)
# Assuming a value of at least 1 for each.
vp_bonus_humans = 0
if credits > 0: vp_bonus_humans += 1
if productivity > 0: vp_bonus_humans += 1
if discovery > 0: vp_bonus_humans += 1
if influence > 0: vp_bonus_humans += 1
    
# Us'ud: +2 VP per Relic
vp_bonus_usud = relics * 2

# Legarchaea: +1 VP per Colonized Planet
vp_bonus_legarchaea = colonized_planets * 1

# Sum all VP sources
all_vp_sources = {
    "Base Planets": vp_from_planets_base,
    "Technology": vp_from_tech,
    "Military": vp_from_military,
    "Base Relics": vp_from_relics_base,
    "Prestige": vp_from_prestige,
    "Chaeilki Bonus": vp_bonus_chaeilki,
    "Humans Bonus": vp_bonus_humans,
    "Us'ud Bonus": vp_bonus_usud,
    "Legarchaea Bonus": vp_bonus_legarchaea
}

total_vp = sum(all_vp_sources.values())

# Create the final equation string
equation_parts = [str(v) for v in all_vp_sources.values()]
equation_string = " + ".join(equation_parts)

print("Victory Point Calculation:")
print("-" * 25)
for source, value in all_vp_sources.items():
    print(f"{source:<20}: {value} VP")
print("-" * 25)
print(f"Total VP Equation: {equation_string} = {total_vp}")
print(f"Final Score: {total_vp} VP")
