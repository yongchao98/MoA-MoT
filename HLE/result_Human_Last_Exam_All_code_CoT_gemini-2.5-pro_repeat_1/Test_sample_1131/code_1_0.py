# Step 1: Define initial values from the game state
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

# The problem states three technologies were researched
technologies_researched = 3

# Step 2: Calculate VPs from various sources based on game rules

# Total number of planets
total_planets = colonized_planets + developed_planets + conquered_planets

# VP from Chaeilki Alliance (1 VP per planet)
chaeilki_vp = total_planets * 1

# VP from Humans Alliance (1 VP per 2 technologies)
humans_vp = technologies_researched // 2

# VP from Us'ud Alliance (1 VP per 2 cruisers)
usud_vp = cruisers // 2

# VP from Legarchaea Ideology (1 VP per 3 planets)
legarchaea_vp = total_planets // 3

# VP from Prestige (1 VP per prestige point)
prestige_vp = prestige * 1

# VP from Relics (assuming 3 VP per relic, a standard rule)
relic_vp = relics * 3

# Bonus VP from Developed Planets (1 bonus VP per developed planet)
developed_planet_bonus_vp = developed_planets * 1

# Step 3: Sum the VPs for the total score
total_vp = chaeilki_vp + humans_vp + usud_vp + legarchaea_vp + prestige_vp + relic_vp + developed_planet_bonus_vp

# Step 4: Print the final equation and result
print(f"VP from Chaeilki (planets): {chaeilki_vp} + VP from Humans (techs): {humans_vp} + VP from Us'ud (cruisers): {usud_vp} + VP from Legarchaea (planets): {legarchaea_vp} + VP from Prestige: {prestige_vp} + VP from Relics: {relic_vp} + VP from Developed Planet Bonus: {developed_planet_bonus_vp} = {total_vp}")
