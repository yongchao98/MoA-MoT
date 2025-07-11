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
researched_techs_count = 3 # Terraforming, Advanced Terraforming, Planetary Shields

# --- VP Calculation ---

# 1. Base VP from resources, planets, and research
vp_prestige = prestige * 1
vp_relics_base = relics * 1
vp_cruisers = cruisers * 1
vp_colonized_base = colonized_planets * 2
vp_developed_base = developed_planets * 3
vp_conquered = conquered_planets * 1
vp_research = researched_techs_count * 1

# 2. VP from Legarchaea ideology
# 1 VP for every 3 combined credits, productivity, discovery, and influence
total_resources = credits + productivity + discovery + influence
vp_ideology = total_resources // 3

# 3. VP from Alliances
# Chaeilki: 2 VP for each relic
vp_alliance_chaeilki = relics * 2
# Humans: 1 VP for each developed planet
vp_alliance_humans = developed_planets * 1
# Us'ud: 2 VP for each colonized planet
vp_alliance_usud = colonized_planets * 2

# 4. Summing it all up
total_vp = (vp_prestige + vp_relics_base + vp_cruisers +
             vp_colonized_base + vp_developed_base + vp_conquered +
             vp_research + vp_ideology + vp_alliance_chaeilki +
             vp_alliance_humans + vp_alliance_usud)

# --- Output the results ---
print("Calculating Victory Points Breakdown:")
print(f"Prestige: {prestige} * 1 VP/prestige = {vp_prestige} VP")
print(f"Relics (base): {relics} * 1 VP/relic = {vp_relics_base} VP")
print(f"Cruisers: {cruisers} * 1 VP/cruiser = {vp_cruisers} VP")
print(f"Colonized Planets (base): {colonized_planets} * 2 VP/planet = {vp_colonized_base} VP")
print(f"Developed Planets (base): {developed_planets} * 3 VP/planet = {vp_developed_base} VP")
print(f"Conquered Planets: {conquered_planets} * 1 VP/planet = {vp_conquered} VP")
print(f"Researched Techs: {researched_techs_count} * 1 VP/tech = {vp_research} VP")
print(f"Ideology (Legarchaea): ({credits}+{productivity}+{discovery}+{influence}) / 3 = {total_resources} / 3 = {vp_ideology} VP")
print(f"Alliance (Chaeilki): {relics} relics * 2 VP/relic = {vp_alliance_chaeilki} VP")
print(f"Alliance (Humans): {developed_planets} developed planet * 1 VP/planet = {vp_alliance_humans} VP")
print(f"Alliance (Us'ud): {colonized_planets} colonized planets * 2 VP/planet = {vp_alliance_usud} VP")
print("\n---")
print("Final VP Calculation:")
# In the final print, we still need to output each number in the equation.
print(f"Total VP = {vp_prestige} + {vp_relics_base} + {vp_cruisers} + {vp_colonized_base} + {vp_developed_base} + {vp_conquered} + {vp_research} + {vp_ideology} + {vp_alliance_chaeilki} + {vp_alliance_humans} + {vp_alliance_usud} = {total_vp}")
<<<42>>>