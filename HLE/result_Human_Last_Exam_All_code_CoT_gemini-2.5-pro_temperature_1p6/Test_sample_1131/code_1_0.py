# Game State Variables
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
# Techs are flags
terraforming_researched = True
advanced_terraforming_researched = True
planetary_shields_researched = True # This tech does not provide VP

# --- VP Calculation ---

# 1. Alliance VP
chaeilki_vp = relics * 2  # 2 VP per relic for Chaeilki
humans_vp = developed_planets * 1  # 1 VP per developed planet for Humans
usud_vp = conquered_planets * 1  # 1 VP per conquered planet for Us'ud

# 2. Ideology VP (Legarchaea)
legarchaea_vp = 0
if prestige >= 7:
    legarchaea_vp = 6
elif prestige >= 5:
    legarchaea_vp = 3
elif prestige >= 3:
    legarchaea_vp = 1

# 3. Standard VP
colonized_planets_vp = colonized_planets * 1
developed_planets_vp = developed_planets * 2
conquered_planets_vp = conquered_planets * 2
cruisers_vp = cruisers // 2  # 1 VP per 2 cruisers
relics_vp = relics * 1 # Standard 1 VP per relic
tech_vp = 0
if terraforming_researched:
    tech_vp += 1
if advanced_terraforming_researched:
    tech_vp += 2

# 4. Resource VP
total_resources = credits + productivity + discovery
resource_vp = total_resources // 10

# --- Total VP Calculation ---
total_vp = (
    chaeilki_vp + humans_vp + usud_vp + legarchaea_vp +
    colonized_planets_vp + developed_planets_vp + conquered_planets_vp +
    cruisers_vp + relics_vp + tech_vp + resource_vp
)

# --- Output ---
print("Victory Point Calculation Breakdown:")
print(f"- Alliance (Chaeilki): {relics} relics * 2 = {chaeilki_vp} VP")
print(f"- Alliance (Humans): {developed_planets} developed planets * 1 = {humans_vp} VP")
print(f"- Alliance (Us'ud): {conquered_planets} conquered planet * 1 = {usud_vp} VP")
print(f"- Ideology (Legarchaea): {prestige} prestige = {legarchaea_vp} VP")
print(f"- Standard (Colonized Planets): {colonized_planets} planets * 1 = {colonized_planets_vp} VP")
print(f"- Standard (Developed Planets): {developed_planets} planet * 2 = {developed_planets_vp} VP")
print(f"- Standard (Conquered Planets): {conquered_planets} planet * 2 = {conquered_planets_vp} VP")
print(f"- Standard (Cruisers): {cruisers} cruisers // 2 = {cruisers_vp} VP")
print(f"- Standard (Relics): {relics} relics * 1 = {relics_vp} VP")
print(f"- Standard (Techs): Terraforming + Advanced Terraforming = {1 if terraforming_researched else 0} + {2 if advanced_terraforming_researched else 0} = {tech_vp} VP")
print(f"- Standard (Resources): ({credits} + {productivity} + {discovery}) // 10 = {resource_vp} VP")
print("\nFinal VP Equation:")
print(f"{chaeilki_vp} + {humans_vp} + {usud_vp} + {legarchaea_vp} + {colonized_planets_vp} + {developed_planets_vp} + {conquered_planets_vp} + {cruisers_vp} + {relics_vp} + {tech_vp} + {resource_vp} = {total_vp}")
print(f"\nTotal Victory Points: {total_vp}")
