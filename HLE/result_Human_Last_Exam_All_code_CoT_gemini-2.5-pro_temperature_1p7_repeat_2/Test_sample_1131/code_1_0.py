import math

# Player's assets
developed_planets = 1
colonized_planets = 3
conquered_planets = 1
cruisers = 5
prestige = 5
relics = 2

# VP calculation multipliers
chaeilki_vp_per_developed = 2
human_vp_per_colonized = 1
usud_vp_per_conquered = 1
prestige_vp_multiplier = 1
relic_vp_multiplier = 2

# Calculate VP from each source
chaeilki_vp = developed_planets * chaeilki_vp_per_developed
human_vp = colonized_planets * human_vp_per_colonized
usud_vp = conquered_planets * usud_vp_per_conquered
legarchaea_vp = cruisers // 2  # Integer division rounds down
prestige_vp = prestige * prestige_vp_multiplier
relics_vp = relics * relic_vp_multiplier

# Calculate total VP
total_vp = chaeilki_vp + human_vp + usud_vp + legarchaea_vp + prestige_vp + relics_vp

# Print the breakdown
print("Calculating Total Victory Points:")
print(f"Chaeilki Faction VP ({developed_planets} developed planet(s) * {chaeilki_vp_per_developed} VP): {chaeilki_vp} VP")
print(f"Human Faction VP ({colonized_planets} colonized planet(s) * {human_vp_per_colonized} VP): {human_vp} VP")
print(f"Us'ud Faction VP ({conquered_planets} conquered planet(s) * {usud_vp_per_conquered} VP): {usud_vp} VP")
print(f"Legarchaea Ideology VP ({cruisers} cruisers / 2): {legarchaea_vp} VP")
print(f"Prestige VP ({prestige} prestige * {prestige_vp_multiplier} VP): {prestige_vp} VP")
print(f"Relics VP ({relics} relics * {relic_vp_multiplier} VP): {relics_vp} VP")
print("-" * 20)

# Print the final equation and result
print("Final Calculation:")
print(f"{chaeilki_vp} + {human_vp} + {usud_vp} + {legarchaea_vp} + {prestige_vp} + {relics_vp} = {total_vp}")
print(f"Total VP: {total_vp}")
