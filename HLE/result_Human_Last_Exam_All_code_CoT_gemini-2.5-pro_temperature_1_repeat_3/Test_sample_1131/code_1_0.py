# Player stats
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
researched_techs = 3 # Terraforming, Advanced Terraforming, Planetary Shields

# --- VP Calculation ---

# General VP
credits_vp = credits // 5
productivity_vp = productivity // 3
discovery_vp = discovery // 3
influence_vp = influence // 3
prestige_vp = prestige
relics_vp = relics * 3
cruisers_vp = cruisers
total_planets = colonized_planets + developed_planets + conquered_planets
planets_vp = total_planets * 2
developed_planet_bonus_vp = developed_planets * 1
tech_vp = researched_techs * 1

# Alliance VP
chaeilki_vp = relics * 2
humans_vp = total_planets
usud_vp = productivity // 3

# Ideology VP
legarchaea_vp = developed_planets * 1

# Total VP
total_vp = (credits_vp + productivity_vp + discovery_vp + influence_vp +
            prestige_vp + relics_vp + cruisers_vp + planets_vp +
            developed_planet_bonus_vp + tech_vp + chaeilki_vp +
            humans_vp + usud_vp + legarchaea_vp)

# --- Output the results ---

print("VP Calculation Breakdown:")
print(f"Credits: {credits_vp} VP")
print(f"Productivity: {productivity_vp} VP")
print(f"Discovery: {discovery_vp} VP")
print(f"Influence: {influence_vp} VP")
print(f"Prestige: {prestige_vp} VP")
print(f"Relics: {relics_vp} VP")
print(f"Cruisers: {cruisers_vp} VP")
print(f"Planets (Colonized/Developed/Conquered): {planets_vp} VP")
print(f"Developed Planet Bonus: {developed_planet_bonus_vp} VP")
print(f"Researched Technologies: {tech_vp} VP")
print("-" * 20)
print("Alliance VP:")
print(f"Chaeilki (Relics): {chaeilki_vp} VP")
print(f"Humans (Planets): {humans_vp} VP")
print(f"Us'ud (Productivity): {usud_vp} VP")
print("-" * 20)
print("Ideology VP:")
print(f"Legarchaea (Developed Planets): {legarchaea_vp} VP")
print("=" * 20)
print("Final Score Equation:")
print(f"{credits_vp} (Credits) + {productivity_vp} (Productivity) + {discovery_vp} (Discovery) + {influence_vp} (Influence) + "
      f"{prestige_vp} (Prestige) + {relics_vp} (Relics) + {cruisers_vp} (Cruisers) + "
      f"{planets_vp} (Planets) + {developed_planet_bonus_vp} (Dev. Bonus) + {tech_vp} (Tech) + "
      f"{chaeilki_vp} (Chaeilki) + {humans_vp} (Humans) + {usud_vp} (Us'ud) + "
      f"{legarchaea_vp} (Legarchaea) = {total_vp}")
print(f"\nTotal VP: {total_vp}")
