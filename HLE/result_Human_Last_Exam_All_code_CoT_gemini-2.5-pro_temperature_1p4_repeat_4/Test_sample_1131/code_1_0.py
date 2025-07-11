# Game state variables
num_races_in_alliance = 3
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

# VP values for different items
# Alliances: 1 race = 1 VP, 2 races = 3 VP, 3 races = 6 VP, 4 races = 10 VP
if num_races_in_alliance == 1:
    alliance_vp = 1
elif num_races_in_alliance == 2:
    alliance_vp = 3
elif num_races_in_alliance == 3:
    alliance_vp = 6
elif num_races_in_alliance == 4:
    alliance_vp = 10
else:
    alliance_vp = 0

# Ideology (Legarchaea): 2 VP per relic
ideology_vp = relics * 2

# Resource sets: 1 VP per set of (credit, productivity, discovery, influence)
resource_set_vp = min(credits, productivity, discovery, influence) * 1

# Prestige: 1 VP per prestige point
prestige_vp = prestige * 1

# Fleet: 1 VP per cruiser
fleet_vp = cruisers * 1

# Planets: Colonized=1 VP, Developed=2 VP, Conquered=1 VP
colonized_planet_vp = colonized_planets * 1
developed_planet_vp = developed_planets * 2
conquered_planet_vp = conquered_planets * 1

# Technologies: Terraforming=1, Advanced Terraforming=2, Planetary Shields=1
tech_terraforming_vp = 1
tech_advanced_terraforming_vp = 2
tech_planetary_shields_vp = 1

# Calculate total VP
total_vp = (alliance_vp + ideology_vp + resource_set_vp + prestige_vp + fleet_vp +
            colonized_planet_vp + developed_planet_vp + conquered_planet_vp +
            tech_terraforming_vp + tech_advanced_terraforming_vp + tech_planetary_shields_vp)

# Print the detailed equation
print(f"Calculating total Victory Points:")
print(f"{alliance_vp} (Alliance) + {ideology_vp} (Ideology) + {resource_set_vp} (Resources) + {prestige_vp} (Prestige) + {fleet_vp} (Fleet) + {colonized_planet_vp} (Colonized) + {developed_planet_vp} (Developed) + {conquered_planet_vp} (Conquered) + {tech_terraforming_vp} (Terraforming) + {tech_advanced_terraforming_vp} (Adv. Terraforming) + {tech_planetary_shields_vp} (Planetary Shields)")
print(f"= {total_vp} VP")
