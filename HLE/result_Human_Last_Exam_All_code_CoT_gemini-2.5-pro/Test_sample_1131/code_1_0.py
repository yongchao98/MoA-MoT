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

# VP values based on Age of Galaxy rules
vp_per_relic = 2
vp_per_cruiser = 1
vp_per_colonized = 1
vp_per_developed = 2
vp_per_conquered = 1

# Alliance VPs
chaeilki_vp = 4
humans_vp = 2
usud_vp = 2
alliance_vp = chaeilki_vp + humans_vp + usud_vp
print(f"Alliance VP (Chaeilki + Humans + Us'ud): {chaeilki_vp} + {humans_vp} + {usud_vp} = {alliance_vp}")

# Ideology VP
legarchaea_vp = 4
print(f"Ideology VP (Legarchaea): {legarchaea_vp}")

# Resource VPs (Credits, Productivity, Discovery, Influence, Prestige)
total_resources = credits + productivity + discovery + influence + prestige
resource_vp = total_resources // 5
print(f"Resource VP (({credits} + {productivity} + {discovery} + {influence} + {prestige}) / 5): {resource_vp}")

# Relic VPs
relic_vp = relics * vp_per_relic
print(f"Relic VP ({relics} * {vp_per_relic}): {relic_vp}")

# Military VP
military_vp = cruisers * vp_per_cruiser
print(f"Military VP ({cruisers} * {vp_per_cruiser}): {military_vp}")

# Planet VP
planet_vp = (colonized_planets * vp_per_colonized) + \
            (developed_planets * vp_per_developed) + \
            (conquered_planets * vp_per_conquered)
print(f"Planet VP (({colonized_planets} * {vp_per_colonized}) + ({developed_planets} * {vp_per_developed}) + ({conquered_planets} * {vp_per_conquered})): {planet_vp}")

# Technology VP
# Terraforming, Advanced Terraforming, and Planetary Shields do not provide direct VP.
tech_vp = 0
print(f"Technology VP: {tech_vp}")

# Total VP
total_vp = alliance_vp + legarchaea_vp + resource_vp + relic_vp + military_vp + planet_vp + tech_vp
print("\n--- Total VP Calculation ---")
print(f"Final Score = {alliance_vp} (Alliances) + {legarchaea_vp} (Ideology) + {resource_vp} (Resources) + {relic_vp} (Relics) + {military_vp} (Military) + {planet_vp} (Planets) + {tech_vp} (Tech)")
print(f"Total VP = {total_vp}")

print("<<<31>>>")