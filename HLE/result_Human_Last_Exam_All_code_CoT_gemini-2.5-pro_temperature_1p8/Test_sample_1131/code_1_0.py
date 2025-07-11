import math

# Player's assets
prestige = 5
colonized_planets = 3
relics = 2
cruisers = 5
developed_planets = 1
conquered_planets = 1

# VP calculation
total_vp = 0

# 1. VP from Alliances
# Humans: 1 VP for every 2 Prestige
vp_humans = prestige // 2
# Chaeilki: 2 VP for every Colonized Planet
vp_chaeilki = colonized_planets * 2
# Us'ud: 1 VP for every Relic
vp_usud = relics * 1

# 2. VP from Ideology
# Legarchaea: 1 VP for every 2 Cruisers
vp_legarchaea = cruisers // 2

# 3. VP from Planets
# Colonized Planet: 1 VP each
vp_colonized = colonized_planets * 1
# Developed Planet: 2 VP each
vp_developed = developed_planets * 2
# Conquered Planet: 1 VP each
vp_conquered = conquered_planets * 1

# 4. VP from Research
# Planetary Shields: 1 VP per developed planet
vp_planetary_shields = developed_planets * 1

# Summing it all up
total_vp = (vp_humans + vp_chaeilki + vp_usud + vp_legarchaea + 
            vp_colonized + vp_developed + vp_conquered + vp_planetary_shields)

# Print the final equation
print("Calculating Total Victory Points:")
print(f"VP from Humans (1 per 2 Prestige): {vp_humans}")
print(f"VP from Chaeilki (2 per Colonized Planet): {vp_chaeilki}")
print(f"VP from Us'ud (1 per Relic): {vp_usud}")
print(f"VP from Legarchaea (1 per 2 Cruisers): {vp_legarchaea}")
print(f"VP from Colonized Planets (1 each): {vp_colonized}")
print(f"VP from Developed Planets (2 each): {vp_developed}")
print(f"VP from Conquered Planets (1 each): {vp_conquered}")
print(f"VP from Planetary Shields (1 per Developed Planet): {vp_planetary_shields}")
print("-" * 20)
print(f"Final Equation: {vp_humans} + {vp_chaeilki} + {vp_usud} + {vp_legarchaea} + {vp_colonized} + {vp_developed} + {vp_conquered} + {vp_planetary_shields} = {total_vp}")
print(f"Total VP: {total_vp}")

<<<19>>>