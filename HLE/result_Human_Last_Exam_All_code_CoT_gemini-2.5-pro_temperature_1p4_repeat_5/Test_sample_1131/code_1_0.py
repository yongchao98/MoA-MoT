import sys

# Player's assets
influence = 3
prestige = 5
relics = 2
cruisers = 5
colonized_planets = 3
developed_planets = 1
conquered_planets = 1

# --- VP Calculation ---

# 1. Alliance VPs
# Chaeilki: 3 VP for each planet with 1 or more Cruisers.
# We have 5 total planets and 5 cruisers, so we can place one on each.
total_planets = colonized_planets + developed_planets + conquered_planets
planets_with_cruisers = min(total_planets, cruisers)
chaeilki_vp = planets_with_cruisers * 3

# Humans: 2 VP for each Developed planet, 1 VP for each Colonized or Conquered planet.
human_vp = (developed_planets * 2) + ((colonized_planets + conquered_planets) * 1)

# Us'ud: 1 VP for each Cruiser.
usud_vp = cruisers * 1

# 2. Ideology VPs
# Legarchaea (Major): 1 VP per Relic and per Influence point.
legarchaea_vp = relics + influence

# 3. General VPs
# Prestige: 1 VP per Prestige point.
prestige_vp = prestige * 1

# Developed Planets: 1 VP for each developed planet (general rule).
developed_planet_general_vp = developed_planets * 1

# Total VP Calculation
total_vp = chaeilki_vp + human_vp + usud_vp + legarchaea_vp + prestige_vp + developed_planet_general_vp

# --- Output ---
print("VP Calculation Breakdown:")
print(f"Chaeilki Alliance ({planets_with_cruisers} planets with cruisers * 3 VP): {chaeilki_vp} VP")
print(f"Human Alliance ({developed_planets} dev * 2 + {colonized_planets + conquered_planets} other * 1): {human_vp} VP")
print(f"Us'ud Alliance ({cruisers} cruisers * 1 VP): {usud_vp} VP")
print(f"Legarchaea Ideology ({relics} relics + {influence} influence): {legarchaea_vp} VP")
print(f"Prestige ({prestige} prestige * 1 VP): {prestige_vp} VP")
print(f"General Rule ({developed_planets} developed planet * 1 VP): {developed_planet_general_vp} VP")
print("-" * 20)
print("Final Equation:")
print(f"{chaeilki_vp} + {human_vp} + {usud_vp} + {legarchaea_vp} + {prestige_vp} + {developed_planet_general_vp} = {total_vp}")
print("-" * 20)
print(f"Total VP: {total_vp}")

# This line is for the final answer extraction.
sys.stdout.write(f'<<<{total_vp}>>>')