# Player's assets
num_relics = 2
num_colonized_planets = 3
num_developed_planets = 1
num_conquered_planets = 1
prestige = 5

# VP values
vp_per_colonized = 2
vp_per_developed = 3
vp_per_conquered = 1

# VP from planets (base values)
colonized_planet_vp = num_colonized_planets * vp_per_colonized
developed_planet_vp = num_developed_planets * vp_per_developed
conquered_planet_vp = num_conquered_planets * vp_per_conquered

# VP from Alliance bonuses
# Chaeilki: +1 VP per Relic
chaeilki_bonus_vp = num_relics * 1
# Humans: +1 VP per Colonized Planet
humans_bonus_vp = num_colonized_planets * 1
# Us'ud: +1 VP per Conquered Planet
usud_bonus_vp = num_conquered_planets * 1

# VP from Ideology bonus
# Legarchaea: +2 VP per Developed Planet
legarchaea_bonus_vp = num_developed_planets * 2

# VP from Prestige
prestige_vp = prestige // 5

# Calculate total VP
total_vp = (
    colonized_planet_vp +
    developed_planet_vp +
    conquered_planet_vp +
    chaeilki_bonus_vp +
    humans_bonus_vp +
    usud_bonus_vp +
    legarchaea_bonus_vp +
    prestige_vp
)

# Print the breakdown of the final score
print("Calculating total Victory Points (VP):")
print(f"- From {num_colonized_planets} Colonized Planets: {colonized_planet_vp} VP")
print(f"- From {num_developed_planets} Developed Planet: {developed_planet_vp} VP")
print(f"- From {num_conquered_planets} Conquered Planet: {conquered_planet_vp} VP")
print(f"- From Chaeilki Alliance ({num_relics} Relics): {chaeilki_bonus_vp} VP")
print(f"- From Human Alliance ({num_colonized_planets} Colonized Planets): {humans_bonus_vp} VP")
print(f"- From Us'ud Alliance ({num_conquered_planets} Conquered Planet): {usud_bonus_vp} VP")
print(f"- From Legarchaea Ideology ({num_developed_planets} Developed Planet): {legarchaea_bonus_vp} VP")
print(f"- From {prestige} Prestige: {prestige_vp} VP")
print("\nFinal VP Calculation:")
print(f"{colonized_planet_vp} (Colonized) + {developed_planet_vp} (Developed) + {conquered_planet_vp} (Conquered) + {chaeilki_bonus_vp} (Chaeilki) + {humans_bonus_vp} (Humans) + {usud_bonus_vp} (Us'ud) + {legarchaea_bonus_vp} (Legarchaea) + {prestige_vp} (Prestige) = {total_vp}")
print(f"\nTotal VP: {total_vp}")

<<<19>>>