import collections

# --- Player's Assets ---
# Resources
prestige = 5
relics = 2
cruisers = 5
# Planets
colonized_planets = 3
developed_planets = 1
conquered_planets = 1
# Techs (with their VP values)
researched_techs = {
    "Terraforming": 2,
    "Advanced Terraforming": 3,
    "Planetary Shields": 1,
}

# --- VP Calculation ---
total_vp = 0
# Using an ordered dictionary to maintain the order of calculation for printing
equation_parts = collections.OrderedDict()

# 1. Base Planet VPs
vp_from_colonized = colonized_planets * 2
equation_parts["Colonized Planets"] = (f"{colonized_planets} * 2", vp_from_colonized)

vp_from_developed = developed_planets * 3
equation_parts["Developed Planets"] = (f"{developed_planets} * 3", vp_from_developed)

vp_from_conquered = conquered_planets * 1
equation_parts["Conquered Planets"] = (f"{conquered_planets} * 1", vp_from_conquered)

# 2. Researched Technology VPs
vp_from_tech = sum(researched_techs.values())
tech_vp_strings = [str(vp) for vp in researched_techs.values()]
equation_parts["Technologies"] = (f"({' + '.join(tech_vp_strings)})", vp_from_tech)

# 3. Prestige VPs
vp_from_prestige = prestige * 1
equation_parts["Prestige"] = (f"{prestige} * 1", vp_from_prestige)

# 4. Base Relic VPs
vp_from_relics_base = relics * 3
equation_parts["Relics (Base)"] = (f"{relics} * 3", vp_from_relics_base)

# 5. Faction-Specific VPs
# Chaeilki: 2 VP per Relic
vp_from_chaeilki = relics * 2
equation_parts["Chaeilki (Relics)"] = (f"{relics} * 2", vp_from_chaeilki)

# Humans: 3 VP per Developed Planet
vp_from_humans = developed_planets * 3
equation_parts["Humans (Developed Planets)"] = (f"{developed_planets} * 3", vp_from_humans)

# Us'ud: 1 VP per Cruiser
vp_from_usud = cruisers * 1
equation_parts["Us'ud (Cruisers)"] = (f"{cruisers} * 1", vp_from_usud)

# 6. Ideology-Specific VPs
# Legarchaea: 2 VP per Colonized or Conquered Planet
planets_for_ideology = colonized_planets + conquered_planets
vp_from_legarchaea = planets_for_ideology * 2
equation_parts["Legarchaea (Planets)"] = (f"({colonized_planets} + {conquered_planets}) * 2", vp_from_legarchaea)


# --- Final Output ---
print("Calculating total Victory Points (VP)...")
print("-" * 30)
# Build the strings for the final equation and sum the total VP
final_equation_str_parts = []
for source, (calc_str, vp) in equation_parts.items():
    print(f"VP from {source}: {vp}")
    final_equation_str_parts.append(calc_str)
    total_vp += vp
print("-" * 30)

# Print the final combined equation
final_equation_str = " + ".join(final_equation_str_parts)
print(f"Total VP = {final_equation_str}")
print(f"Total VP = {total_vp}")
print("<<<47>>>")