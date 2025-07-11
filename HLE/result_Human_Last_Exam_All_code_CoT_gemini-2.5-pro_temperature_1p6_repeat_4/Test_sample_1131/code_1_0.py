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
# Researched Techs
terraforming_tech = 1
advanced_terraforming_tech = 1
planetary_shields_tech = 1

# --- VP Calculation ---
total_vp = 0
calculation_steps = []

# 1. Alliance VPs
chaeilki_vp = relics * 2
total_vp += chaeilki_vp
calculation_steps.append(str(chaeilki_vp))

humans_vp = colonized_planets * 2
total_vp += humans_vp
calculation_steps.append(str(humans_vp))

usud_vp = developed_planets * 5
total_vp += usud_vp
calculation_steps.append(str(usud_vp))

# 2. Major Ideology VPs
total_planets = colonized_planets + developed_planets + conquered_planets
legarchaea_vp = total_planets * 1
total_vp += legarchaea_vp
calculation_steps.append(str(legarchaea_vp))

# 3. Standard VPs
colonized_vp_base = colonized_planets * 1
total_vp += colonized_vp_base
calculation_steps.append(str(colonized_vp_base))

developed_vp_base = developed_planets * 3
total_vp += developed_vp_base
calculation_steps.append(str(developed_vp_base))

conquered_vp_base = conquered_planets * 2
total_vp += conquered_vp_base
calculation_steps.append(str(conquered_vp_base))

credit_vp = (credits // 7) * 1
total_vp += credit_vp
calculation_steps.append(str(credit_vp))

influence_vp = (influence // 3) * 1
total_vp += influence_vp
calculation_steps.append(str(influence_vp))

prestige_vp = (prestige // 5) * 1
total_vp += prestige_vp
calculation_steps.append(str(prestige_vp))

relics_vp_base = relics * 3
total_vp += relics_vp_base
calculation_steps.append(str(relics_vp_base))

cruisers_vp = cruisers * 1
total_vp += cruisers_vp
calculation_steps.append(str(cruisers_vp))

# 4. Technology VPs
terraforming_vp = 1 * terraforming_tech
total_vp += terraforming_vp
calculation_steps.append(str(terraforming_vp))

advanced_terraforming_vp = 2 * advanced_terraforming_tech
total_vp += advanced_terraforming_vp
calculation_steps.append(str(advanced_terraforming_vp))

planetary_shields_vp = 2 * planetary_shields_tech
total_vp += planetary_shields_vp
calculation_steps.append(str(planetary_shields_vp))

# --- Output the result ---
print("Calculating Total Victory Points (VP):")
print(f"Alliance - Chaeilki (2 relics * 2 VP): {chaeilki_vp}")
print(f"Alliance - Humans (3 colonized planets * 2 VP): {humans_vp}")
print(f"Alliance - Us'ud (1 developed planet * 5 VP): {usud_vp}")
print(f"Ideology - Legarchaea (5 controlled planets * 1 VP): {legarchaea_vp}")
print(f"Planets - Colonized (3 planets * 1 VP): {colonized_vp_base}")
print(f"Planets - Developed (1 planet * 3 VP): {developed_vp_base}")
print(f"Planets - Conquered (1 planet * 2 VP): {conquered_vp_base}")
print(f"Resources - Credits (1 set of 7): {credit_vp}")
print(f"Resources - Influence (1 set of 3): {influence_vp}")
print(f"Resources - Prestige (1 set of 5): {prestige_vp}")
print(f"Items - Relics (2 relics * 3 VP): {relics_vp_base}")
print(f"Units - Cruisers (5 cruisers * 1 VP): {cruisers_vp}")
print(f"Tech - Terraforming: {terraforming_vp}")
print(f"Tech - Advanced Terraforming: {advanced_terraforming_vp}")
print(f"Tech - Planetary Shields: {planetary_shields_vp}")
print("-" * 20)
final_equation = " + ".join(calculation_steps)
print(f"Total VP = {final_equation} = {total_vp}")
print("<<<47>>>")