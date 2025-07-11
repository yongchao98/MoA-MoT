# Player's assets
alliances_count = 3
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

# VP Calculations

# 1. Ideology (Legarchaea) VP: 2 VP per alliance
ideology_vp = alliances_count * 2

# 2. Alliance-Specific VP
# Chaeilki: 1 VP per colonized planet
chaeilki_vp = colonized_planets * 1
# Humans: 1 VP per prestige
humans_vp = prestige * 1
# Us'ud: 2 VP per relic
usud_vp = relics * 2

# 3. Standard VP from Assets
# Resource sets: 1 VP per set of (credit, productivity, discovery, influence)
resource_sets = min(credits, productivity, discovery, influence)
resource_vp = resource_sets * 1
# Standard Prestige VP: 1 VP per 2 prestige
prestige_vp = prestige // 2
# Standard Cruiser VP: 1 VP per 2 cruisers
cruiser_vp = cruisers // 2
# Standard Relic VP: 1 VP per relic
relic_vp = relics * 1

# 4. Standard VP from Planets
# Colonized planets: 1 VP per planet
colonized_vp = colonized_planets * 1
# Developed planets: 2 VP per planet
developed_vp = developed_planets * 2
# Conquered planets: 1 VP per planet
conquered_vp = conquered_planets * 1

# 5. VP from Research
# Terraforming: 1 VP, Advanced Terraforming: 2 VP, Planetary Shields: 1 VP
research_vp = 1 + 2 + 1

# 6. Total VP
total_vp = (ideology_vp + chaeilki_vp + humans_vp + usud_vp +
            resource_vp + prestige_vp + cruiser_vp + relic_vp +
            colonized_vp + developed_vp + conquered_vp + research_vp)

# Print the detailed breakdown and final calculation
print("VP Calculation Breakdown:")
print(f"Ideology (Legarchaea): {ideology_vp} VP")
print(f"Alliance - Chaeilki (from colonized planets): {chaeilki_vp} VP")
print(f"Alliance - Humans (from prestige): {humans_vp} VP")
print(f"Alliance - Us'ud (from relics): {usud_vp} VP")
print(f"Standard - Resource Sets: {resource_vp} VP")
print(f"Standard - Prestige: {prestige_vp} VP")
print(f"Standard - Cruisers: {cruiser_vp} VP")
print(f"Standard - Relics: {relic_vp} VP")
print(f"Standard - Colonized Planets: {colonized_vp} VP")
print(f"Standard - Developed Planets: {developed_vp} VP")
print(f"Standard - Conquered Planets: {conquered_vp} VP")
print(f"Research (Terraforming, Adv. Terraforming, Planetary Shields): {research_vp} VP")
print("-" * 25)
print("Total VP Equation:")
print(f"{ideology_vp} + {chaeilki_vp} + {humans_vp} + {usud_vp} + {resource_vp} + {prestige_vp} + {cruiser_vp} + {relic_vp} + {colonized_vp} + {developed_vp} + {conquered_vp} + {research_vp} = {total_vp}")
print("-" * 25)
print(f"Final Total VP: {total_vp}")