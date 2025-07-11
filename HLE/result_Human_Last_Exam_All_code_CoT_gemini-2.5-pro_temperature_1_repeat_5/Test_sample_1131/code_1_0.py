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
# Technologies are flags, 1 for researched, 0 for not
tech_terraforming = 1
tech_advanced_terraforming = 2 # This tech is worth 2 VP
tech_planetary_shields = 1

total_vp = 0
equation_parts = []

# 1. Alliance VPs
# Chaeilki: 1 VP per Relic
chaeilki_vp = relics * 1
total_vp += chaeilki_vp
equation_parts.append(str(chaeilki_vp))
print(f"Chaeilki VP (from {relics} relics): {chaeilki_vp}")

# Humans: 1 VP per Colonized Planet
humans_vp = colonized_planets * 1
total_vp += humans_vp
equation_parts.append(str(humans_vp))
print(f"Humans VP (from {colonized_planets} colonized planets): {humans_vp}")

# Us'ud: 1 VP per Developed Planet
usud_vp = developed_planets * 1
total_vp += usud_vp
equation_parts.append(str(usud_vp))
print(f"Us'ud VP (from {developed_planets} developed planets): {usud_vp}")

# 2. Ideology VPs
# Legarchaea: 1 VP per Prestige
legarchaea_vp = prestige * 1
total_vp += legarchaea_vp
equation_parts.append(str(legarchaea_vp))
print(f"Legarchaea VP (from {prestige} prestige): {legarchaea_vp}")

# 3. General VPs
# Resources
credits_vp = credits // 7
total_vp += credits_vp
equation_parts.append(str(credits_vp))
print(f"Credits VP (from {credits} credits): {credits_vp}")

productivity_vp = productivity // 5
total_vp += productivity_vp
equation_parts.append(str(productivity_vp))
print(f"Productivity VP (from {productivity} productivity): {productivity_vp}")

discovery_vp = discovery // 4
total_vp += discovery_vp
equation_parts.append(str(discovery_vp))
print(f"Discovery VP (from {discovery} discovery): {discovery_vp}")

influence_vp = influence // 3
total_vp += influence_vp
equation_parts.append(str(influence_vp))
print(f"Influence VP (from {influence} influence): {influence_vp}")

# Cruisers
cruisers_vp = cruisers // 5
total_vp += cruisers_vp
equation_parts.append(str(cruisers_vp))
print(f"Cruisers VP (from {cruisers} cruisers): {cruisers_vp}")

# Planets (General Scoring)
colonized_vp = colonized_planets * 1
total_vp += colonized_vp
equation_parts.append(str(colonized_vp))
print(f"Colonized Planets VP ({colonized_planets} * 1 VP): {colonized_vp}")

developed_vp = developed_planets * 2
total_vp += developed_vp
equation_parts.append(str(developed_vp))
print(f"Developed Planets VP ({developed_planets} * 2 VP): {developed_vp}")

conquered_vp = conquered_planets * 3
total_vp += conquered_vp
equation_parts.append(str(conquered_vp))
print(f"Conquered Planets VP ({conquered_planets} * 3 VP): {conquered_vp}")

# Technologies
total_vp += tech_terraforming
equation_parts.append(str(tech_terraforming))
print(f"Terraforming Tech VP: {tech_terraforming}")

total_vp += tech_advanced_terraforming
equation_parts.append(str(tech_advanced_terraforming))
print(f"Advanced Terraforming Tech VP: {tech_advanced_terraforming}")

total_vp += tech_planetary_shields
equation_parts.append(str(tech_planetary_shields))
print(f"Planetary Shields Tech VP: {tech_planetary_shields}")

# 4. Final Calculation
print("\n--- Total VP Calculation ---")
final_equation = " + ".join(equation_parts)
print(f"Final Equation: {final_equation} = {total_vp}")
print(f"Total VP: {total_vp}")
