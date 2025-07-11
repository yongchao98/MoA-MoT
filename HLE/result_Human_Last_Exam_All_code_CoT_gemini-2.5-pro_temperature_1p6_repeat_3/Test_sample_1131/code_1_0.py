# Player's assets and game state
colonized_planets = 3
developed_planets = 1
conquered_planets = 1
relics = 2
cruisers = 5
prestige = 5

# --- VP Calculations ---

# 1. Ideology: Legarchaea
# 2 VP per Colonized and Developed Planet you control.
vp_legarchaea_planets = 2 * (colonized_planets + developed_planets)
# 3 VP per Conquered Planet you control.
vp_legarchaea_conquered = 3 * conquered_planets

# 2. Alliance Factions
# Chaeilki: 2 VP for each Relic.
vp_chaeilki = 2 * relics
# Humans: 1 VP for each Cruisership.
vp_humans = 1 * cruisers
# Us'ud: 1 VP for each Colonized Planet.
vp_usud = 1 * colonized_planets

# 3. Researched Technologies
# Terraforming: 1 VP for each Colonized Planet.
vp_terraforming = 1 * colonized_planets
# Advanced Terraforming: 1 VP for each Developed Planet.
vp_adv_terraforming = 1 * developed_planets

# 4. Prestige
# Each prestige point is worth 1 VP.
vp_prestige = prestige

# --- Total VP Calculation ---

# Create a list of all VP components to build the equation
vp_components = [
    vp_legarchaea_planets,
    vp_legarchaea_conquered,
    vp_chaeilki,
    vp_humans,
    vp_usud,
    vp_terraforming,
    vp_adv_terraforming,
    vp_prestige,
]

# Calculate the total VP
total_vp = sum(vp_components)

# --- Print the results ---
print("VP Calculation Breakdown:")
print(f"- Legarchaea (Planets): {vp_legarchaea_planets}")
print(f"- Legarchaea (Conquered): {vp_legarchaea_conquered}")
print(f"- Chaeilki (Relics): {vp_chaeilki}")
print(f"- Humans (Cruisers): {vp_humans}")
print(f"- Us'ud (Colonized Planets): {vp_usud}")
print(f"- Terraforming (Colonized Planets): {vp_terraforming}")
print(f"- Advanced Terraforming (Developed Planets): {vp_adv_terraforming}")
print(f"- Prestige: {vp_prestige}")
print("-" * 20)

# Build and print the final equation as requested
equation_string = " + ".join(map(str, vp_components))
print(f"Final VP Equation: {equation_string} = {total_vp}")

print(f"\nTotal Victory Points: {total_vp}")
