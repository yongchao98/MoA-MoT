# Player's assets
num_races_in_alliance = 3
developed_planets = 1
credits = 7
productivity = 5
discovery = 4
influence = 3
prestige = 5
relics = 2
cruisers = 5
colonized_planets = 3
conquered_planets = 1
num_researched_techs = 3

# --- VP Calculation ---

# 1. Alliance VP (3 races = 5 VP)
alliance_vp = 5

# 2. Ideology VP (Legarchaea: 3 VP per developed planet)
ideology_vp = developed_planets * 3

# 3. Resource VP (min of resources * 2)
resource_vp = min(credits, productivity, discovery, influence) * 2

# 4. Prestige VP (1 VP per prestige)
prestige_vp = prestige

# 5. Relic VP (2 VP per relic)
relic_vp = relics * 2

# 6. Fleet VP (1 VP per cruiser)
fleet_vp = cruisers * 1

# 7. Planet VP (Colonized=1, Developed=2, Conquered=1)
planet_vp = (colonized_planets * 1) + (developed_planets * 2) + (conquered_planets * 1)

# 8. Research VP (1 VP per tech)
research_vp = num_researched_techs * 1

# --- Total VP ---
total_vp = alliance_vp + ideology_vp + resource_vp + prestige_vp + relic_vp + fleet_vp + planet_vp + research_vp

# --- Final Output ---
print("The total VP is calculated as follows:")
print(f"{alliance_vp} (Alliance) + {ideology_vp} (Ideology) + {resource_vp} (Resources) + {prestige_vp} (Prestige) + {relic_vp} (Relics) + {fleet_vp} (Fleet) + {planet_vp} (Planets) + {research_vp} (Research) = {total_vp}")