import math

# --- Player's Assets ---
credits = 7
productivity = 5
discovery = 4
influence = 3
prestige = 5
relics = 2
colonized_planets = 3
developed_planets = 1
conquered_planets = 1
researched_techs_count = 3  # Terraforming, Advanced Terraforming, Planetary Shields

# --- VP Calculations (based on common game rules) ---

# Alliances
# Chaeilki: 2 VP per relic
vp_chaeilki = relics * 2
# Humans: 1 VP per 2 researched technologies (rounded down)
vp_humans = researched_techs_count // 2
# Us'ud: 1 VP per 2 planets (rounded down)
total_planets = colonized_planets + developed_planets + conquered_planets
vp_usud = total_planets // 2

# Ideology
# Legarchaea: 5 base VP + 1 VP per developed planet
vp_legarchaea = 5 + (developed_planets * 1)

# Base Planet Score
# 1 VP per colonized/conquered planet, 2 VP per developed planet
vp_planets = (colonized_planets * 1) + (conquered_planets * 1) + (developed_planets * 2)

# Base Resource Score
# Relics: 1 VP each (base value)
vp_relics = relics * 1
# Other Resources: 1 VP per 7 combined credits, productivity, discovery, influence (CPDI)
cpdi_total = credits + productivity + discovery + influence
vp_cpdi = cpdi_total // 7
# Prestige: 1 VP per 3 prestige
vp_prestige = prestige // 3

# --- Final Summation and Output ---

# Total VP
total_vp = vp_chaeilki + vp_humans + vp_usud + vp_legarchaea + vp_planets + vp_relics + vp_cpdi + vp_prestige

# Print the final equation
print("Total VP Calculation:")
print(
    f"{vp_chaeilki} (Chaeilki) + "
    f"{vp_humans} (Humans) + "
    f"{vp_usud} (Us'ud) + "
    f"{vp_legarchaea} (Legarchaea) + "
    f"{vp_planets} (Planets) + "
    f"{vp_relics} (Relics) + "
    f"{vp_cpdi} (CPDI Resources) + "
    f"{vp_prestige} (Prestige) = {total_vp}"
)
print(f"Final score is {total_vp} VP.")

<<<24>>>