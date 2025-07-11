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

# --- VP Calculation ---

# 1. Standard VP
# 1 VP per Prestige
prestige_vp = prestige * 1

# 2 VP per Relic
relics_base_vp = relics * 2

# 1 VP per 3 units of the highest resource (Credits, Productivity, Discovery, or Influence)
# The highest is Credits (7)
resource_vp = credits // 3

# 1 VP per Colonized Planet
colonized_vp = colonized_planets * 1

# 2 VP per Developed Planet
developed_vp = developed_planets * 2

# 2 VP per Conquered Planet
conquered_vp = conquered_planets * 2

# 1 VP per Cruiser
cruiser_vp = cruisers * 1

# 2. Faction VP
# Chaeilki: 1 VP per 2 Productivity
chaeilki_vp = productivity // 2

# Humans: 1 VP per colonized and developed planet
human_vp = colonized_planets + developed_planets

# Us'ud: 2 VP per Relic
usud_vp = relics * 2

# 3. Total VP Calculation
# Sum all VP sources
total_vp = (prestige_vp + relics_base_vp + resource_vp + 
            colonized_vp + developed_vp + conquered_vp + 
            cruiser_vp + chaeilki_vp + human_vp + usud_vp)

# Print the final equation, showing the contribution from each source
print(
    f"{prestige_vp} (Prestige) + "
    f"{relics_base_vp} (Base Relics) + "
    f"{resource_vp} (Credits) + "
    f"{colonized_vp} (Colonized Planets) + "
    f"{developed_vp} (Developed Planets) + "
    f"{conquered_vp} (Conquered Planets) + "
    f"{cruiser_vp} (Cruisers) + "
    f"{chaeilki_vp} (Chaeilki) + "
    f"{human_vp} (Humans) + "
    f"{usud_vp} (Us'ud) = {total_vp}"
)
<<<33>>>