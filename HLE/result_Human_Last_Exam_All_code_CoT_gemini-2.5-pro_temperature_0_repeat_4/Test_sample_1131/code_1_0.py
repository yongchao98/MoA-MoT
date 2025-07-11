import math

# --- Player's Assets ---
# Resources
credits = 7
productivity = 5
discovery = 4
influence = 3
prestige = 5
relics = 2
cruisers = 5

# Planets
colonized_planets = 3
developed_planets = 1
conquered_planets = 1

# Technologies (VP values)
tech_terraforming_vp = 1
tech_advanced_terraforming_vp = 2
tech_planetary_shields_vp = 1

# --- VP Calculation ---
total_vp = 0
calculation_parts = []

# 1. Alliance VPs
chaeilki_vp = cruisers // 2
calculation_parts.append(str(chaeilki_vp))
print(f"Chaeilki Alliance VP (5 cruisers // 2): {chaeilki_vp}")

humans_vp = colonized_planets // 2
calculation_parts.append(str(humans_vp))
print(f"Human Alliance VP (3 colonized planets // 2): {humans_vp}")

usud_vp = relics // 2
calculation_parts.append(str(usud_vp))
print(f"Us'ud Alliance VP (2 relics // 2): {usud_vp}")

# 2. Ideology VP
legarchaea_vp = prestige // 2
calculation_parts.append(str(legarchaea_vp))
print(f"Legarchaea Ideology VP (5 prestige // 2): {legarchaea_vp}")

# 3. Resource VP
resource_vp = min(credits, productivity, discovery, influence)
calculation_parts.append(str(resource_vp))
print(f"Resource VP (min of {credits}, {productivity}, {discovery}, {influence}): {resource_vp}")

# 4. Planet VPs
colonized_planet_vp = colonized_planets * 1
developed_planet_vp = developed_planets * 2
conquered_planet_vp = conquered_planets * 1
total_planet_vp = colonized_planet_vp + developed_planet_vp + conquered_planet_vp
calculation_parts.append(str(total_planet_vp))
print(f"Planet VP ({colonized_planet_vp} from colonized + {developed_planet_vp} from developed + {conquered_planet_vp} from conquered): {total_planet_vp}")

# 5. Technology VPs
total_tech_vp = tech_terraforming_vp + tech_advanced_terraforming_vp + tech_planetary_shields_vp
calculation_parts.append(str(total_tech_vp))
print(f"Technology VP ({tech_terraforming_vp} + {tech_advanced_terraforming_vp} + {tech_planetary_shields_vp}): {total_tech_vp}")

# --- Final Calculation ---
total_vp = chaeilki_vp + humans_vp + usud_vp + legarchaea_vp + resource_vp + total_planet_vp + total_tech_vp
final_equation = " + ".join(calculation_parts)

print("\n--- Total VP ---")
print(f"Final Calculation: {final_equation} = {total_vp}")
print(f"The player's total VP is: {total_vp}")
