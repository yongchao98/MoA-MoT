# Player's assets
credits = 7
productivity = 5
discovery = 4
influence = 3
prestige = 5
relics = 2
cruisers = 5 # Cruisers do not provide VP at the end of the game

colonized_planets = 3
developed_planets = 1
conquered_planets = 1

# VP values from rules
vp_per_colonized = 3
vp_per_developed = 6
vp_per_conquered = 2

vp_terraforming = 1
vp_adv_terraforming = 2
vp_planetary_shields = 1

# --- Calculations ---

# 1. Planet VPs
planet_vp = (colonized_planets * vp_per_colonized) + \
            (developed_planets * vp_per_developed) + \
            (conquered_planets * vp_per_conquered)
print(f"Planets VP: ({colonized_planets} * {vp_per_colonized}) + ({developed_planets} * {vp_per_developed}) + ({conquered_planets} * {vp_per_conquered}) = {planet_vp}")

# 2. Researched Tech VPs
tech_vp = vp_terraforming + vp_adv_terraforming + vp_planetary_shields
print(f"Technology VP: {vp_terraforming} + {vp_adv_terraforming} + {vp_planetary_shields} = {tech_vp}")

# 3. Resource VPs
total_resources = credits + productivity + discovery + influence
resource_vp = total_resources // 5
print(f"Resources VP: ({credits} + {productivity} + {discovery} + {influence}) // 5 = {resource_vp}")

# 4. Prestige VPs
prestige_vp = prestige
print(f"Prestige VP: {prestige} = {prestige_vp}")

# 5. Ideology VP (Legarchaea)
total_planets = colonized_planets + developed_planets + conquered_planets
ideology_vp = total_planets * 1
print(f"Ideology (Legarchaea) VP: {total_planets} planets * 1 = {ideology_vp}")

# 6. Alliance VPs
num_alliances = 3 # Chaeilki, Humans, Us'ud
# Humans: 1 VP per alliance
humans_vp = num_alliances * 1
# Chaeilki: 2 VP per colonized planet
chaeilki_vp = colonized_planets * 2
# Us'ud: 3 VP per relic
usud_vp = relics * 3
alliance_vp = humans_vp + chaeilki_vp + usud_vp
print(f"Alliance VP (Humans + Chaeilki + Us'ud): {humans_vp} + {chaeilki_vp} + {usud_vp} = {alliance_vp}")

# --- Final Score ---
total_vp = planet_vp + tech_vp + resource_vp + prestige_vp + ideology_vp + alliance_vp

print("\n--- Total VP Calculation ---")
print(f"Final Score = {planet_vp} (Planets) + {tech_vp} (Tech) + {resource_vp} (Resources) + {prestige_vp} (Prestige) + {ideology_vp} (Ideology) + {alliance_vp} (Alliances)")
print(f"Total VP = {total_vp}")

print("\n<<<49>>>")