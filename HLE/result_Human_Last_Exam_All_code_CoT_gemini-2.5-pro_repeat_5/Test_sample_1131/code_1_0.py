# Player's assets
prestige = 5
relics = 2
colonized_planets = 3
developed_planets = 1
conquered_planets = 1
# Technologies: Terraforming, Advanced Terraforming, and Planetary Shields
researched_techs = 3

# VP values per item
vp_per_prestige = 1
vp_per_relic = 3
vp_per_colonized_planet = 1
vp_per_conquered_planet = 2

# Ideology (Legarchaea) VP values
vp_per_developed_planet_ideology = 1
vp_per_tech_ideology = 2

# --- VP Calculation ---

# 1. Standard VP
vp_from_prestige = prestige * vp_per_prestige
vp_from_relics = relics * vp_per_relic
vp_from_colonized = colonized_planets * vp_per_colonized_planet
vp_from_conquered = conquered_planets * vp_per_conquered_planet

# 2. Ideology VP (Legarchaea)
vp_from_developed_ideology = developed_planets * vp_per_developed_planet_ideology
vp_from_tech_ideology = researched_techs * vp_per_tech_ideology

# 3. Total VP
total_vp = (
    vp_from_prestige
    + vp_from_relics
    + vp_from_colonized
    + vp_from_conquered
    + vp_from_developed_ideology
    + vp_from_tech_ideology
)

# --- Output ---
print("Calculating Victory Points based on Legarchaea Ideology:")
print(f"Prestige: {prestige} * {vp_per_prestige} VP = {vp_from_prestige} VP")
print(f"Relics: {relics} * {vp_per_relic} VP = {vp_from_relics} VP")
print(f"Colonized Planets: {colonized_planets} * {vp_per_colonized_planet} VP = {vp_from_colonized} VP")
print(f"Conquered Planets: {conquered_planets} * {vp_per_conquered_planet} VP = {vp_from_conquered} VP")
print(f"Developed Planets (Ideology): {developed_planets} * {vp_per_developed_planet_ideology} VP = {vp_from_developed_ideology} VP")
print(f"Researched Techs (Ideology): {researched_techs} * {vp_per_tech_ideology} VP = {vp_from_tech_ideology} VP")
print("-" * 20)
print(f"Total VP = {vp_from_prestige} + {vp_from_relics} + {vp_from_colonized} + {vp_from_conquered} + {vp_from_developed_ideology} + {vp_from_tech_ideology} = {total_vp}")
<<<23>>>