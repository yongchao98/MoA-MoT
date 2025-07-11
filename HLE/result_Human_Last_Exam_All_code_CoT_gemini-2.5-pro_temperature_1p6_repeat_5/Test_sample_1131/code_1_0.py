# Player's game state variables
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

# 1. Ideology VP (Legarchaea)
legarchaea_developed_vp = developed_planets * 2
legarchaea_conquered_vp = conquered_planets * 2
ideology_vp = legarchaea_developed_vp + legarchaea_conquered_vp

# 2. Alliance VP
# Chaeilki (based on colonized planets)
chaeilki_vp = 0
if colonized_planets == 3:
    chaeilki_vp = 3
elif colonized_planets == 4:
    chaeilki_vp = 6
elif colonized_planets == 5:
    chaeilki_vp = 10
elif colonized_planets >= 6:
    chaeilki_vp = 15
elif colonized_planets >= 1:
    chaeilki_vp = 1


# Humans (based on developed planets)
humans_vp = 0
if developed_planets == 1:
    humans_vp = 3
elif developed_planets == 2:
    humans_vp = 7
elif developed_planets >= 3:
    humans_vp = 12

# Us'ud (based on researched technologies)
usud_vp = 0
if researched_techs == 3:
    usud_vp = 2
elif researched_techs == 4:
    usud_vp = 4
elif researched_techs == 5:
    usud_vp = 7
elif researched_techs >= 6:
    usud_vp = 11

alliance_vp = chaeilki_vp + humans_vp + usud_vp

# 3. General VP
resource_sets_vp = min(credits, productivity, discovery, influence)
prestige_vp = prestige // 2
relics_vp = relics * 2
cruisers_vp = cruisers // 2
general_vp = resource_sets_vp + prestige_vp + relics_vp + cruisers_vp

# 4. Total VP
total_vp = ideology_vp + alliance_vp + general_vp

# --- Output the results ---
print("VP Calculation Breakdown:\n")

print(f"1. Ideology VP (Legarchaea):")
print(f"   - VP from {developed_planets} developed planet(s): {legarchaea_developed_vp}")
print(f"   - VP from {conquered_planets} conquered planet(s): {legarchaea_conquered_vp}")
print(f"   => Subtotal: {ideology_vp} VP\n")

print(f"2. Alliance VP:")
print(f"   - Chaeilki ({colonized_planets} colonized planets): {chaeilki_vp} VP")
print(f"   - Humans ({developed_planets} developed planet): {humans_vp} VP")
print(f"   - Us'ud ({researched_techs} researched techs): {usud_vp} VP")
print(f"   => Subtotal: {alliance_vp} VP\n")

print(f"3. General VP:")
print(f"   - Resource Sets (min of {credits}, {productivity}, {discovery}, {influence}): {resource_sets_vp} VP")
print(f"   - Prestige ({prestige} prestige / 2): {prestige_vp} VP")
print(f"   - Relics ({relics} relics * 2): {relics_vp} VP")
print(f"   - Cruisers ({cruisers} cruisers / 2): {cruisers_vp} VP")
print(f"   => Subtotal: {general_vp} VP\n")

print("--- Final Score ---")
print(f"Ideology ({ideology_vp}) + Alliances ({alliance_vp}) + General ({general_vp})")
print(f"Final Equation: {ideology_vp} + {alliance_vp} + {general_vp} = {total_vp}")
print(f"Total VP: {total_vp}")
