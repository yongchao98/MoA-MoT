# Player's assets
relics = 2
prestige = 5
colonized_planets = 3
developed_planets = 1
conquered_planets = 1

# --- VP Calculation ---
total_vp = 0
calculation_steps = []

# 1. Alliance VP
# Chaeilki: 2 VP per Relic
chaeilki_vp = relics * 2
total_vp += chaeilki_vp
calculation_steps.append(str(chaeilki_vp))
print(f"Chaeilki Alliance VP (2 relics * 2 VP): {chaeilki_vp}")

# Humans: 1 VP per Prestige
humans_vp = prestige * 1
total_vp += humans_vp
calculation_steps.append(str(humans_vp))
print(f"Humans Alliance VP (5 prestige * 1 VP): {humans_vp}")

# Us'ud: 1 VP per Colonized Planet
usud_vp = colonized_planets * 1
total_vp += usud_vp
calculation_steps.append(str(usud_vp))
print(f"Us'ud Alliance VP (3 colonized planets * 1 VP): {usud_vp}")

# 2. Ideology VP
# Legarchaea: 1 VP per Relic and 1 VP per Conquered Planet
legarchaea_vp = (relics * 1) + (conquered_planets * 1)
total_vp += legarchaea_vp
calculation_steps.append(str(legarchaea_vp))
print(f"Legarchaea Ideology VP (2 relics * 1 VP + 1 conquered planet * 1 VP): {legarchaea_vp}")

# 3. Planet VP
# 1 VP per Colonized, 2 per Developed, 1 per Conquered
colonized_planet_vp = colonized_planets * 1
developed_planet_vp = developed_planets * 2
conquered_planet_vp = conquered_planets * 1
planet_vp = colonized_planet_vp + developed_planet_vp + conquered_planet_vp
total_vp += planet_vp
calculation_steps.append(str(planet_vp))
print(f"Planet VP (3 colonized * 1 + 1 developed * 2 + 1 conquered * 1): {planet_vp}")

# 4. Technology VP
# Terraforming: 1 VP, Advanced Terraforming: 2 VP, Planetary Shields: 1 VP
tech_vp = 1 + 2 + 1
total_vp += tech_vp
calculation_steps.append(str(tech_vp))
print(f"Technology VP (Terraforming + Adv. Terraforming + Planetary Shields): {tech_vp}")

# 5. Prestige VP
# 1 VP per Prestige
prestige_vp = prestige * 1
total_vp += prestige_vp
calculation_steps.append(str(prestige_vp))
print(f"Prestige VP (5 prestige * 1 VP): {prestige_vp}")

# --- Final Calculation ---
print("\n--- Total VP Calculation ---")
final_equation = " + ".join(calculation_steps)
print(f"Final Equation: {final_equation} = {total_vp}")
print(f"Total VP: {total_vp}")

<<<30>>>