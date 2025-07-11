# Game state variables
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
# Player has 'Terraforming' and 'Advanced Terraforming'
terraforming_techs = 2

# --- VP Calculations ---

# Alliance VPs
chaeilki_vp = terraforming_techs * 1
humans_vp = colonized_planets * 1
usud_vp = relics * 1

# Ideology VPs
legarchaea_vp = developed_planets * 1 + conquered_planets * 1

# General VPs
credits_vp = credits // 3
productivity_vp = productivity // 3
discovery_vp = discovery // 3
influence_vp = influence // 3
prestige_vp = prestige * 1
relics_vp = relics * 1
cruisers_vp = cruisers // 2
planets_vp = (developed_planets * 2) + (conquered_planets * 1)

# --- Total VP ---
total_vp = (
    chaeilki_vp + humans_vp + usud_vp + legarchaea_vp +
    credits_vp + productivity_vp + discovery_vp + influence_vp +
    prestige_vp + relics_vp + cruisers_vp + planets_vp
)

# --- Print the breakdown ---
print("Calculating Total Victory Points:")
print(f"Alliance (Chaeilki): {chaeilki_vp}")
print(f"Alliance (Humans): {humans_vp}")
print(f"Alliance (Us'ud): {usud_vp}")
print(f"Ideology (Legarchaea): {legarchaea_vp}")
print(f"Credits: {credits_vp}")
print(f"Productivity: {productivity_vp}")
print(f"Discovery: {discovery_vp}")
print(f"Influence: {influence_vp}")
print(f"Prestige: {prestige_vp}")
print(f"Relics (General): {relics_vp}")
print(f"Cruisers: {cruisers_vp}")
print(f"Planets (General): {planets_vp}")
print("-" * 20)
print(
    f"Total VP = {chaeilki_vp} + {humans_vp} + {usud_vp} + {legarchaea_vp} + "
    f"{credits_vp} + {productivity_vp} + {discovery_vp} + {influence_vp} + "
    f"{prestige_vp} + {relics_vp} + {cruisers_vp} + {planets_vp} = {total_vp}"
)

print(f"\nThe final score is {total_vp} VP.")
<<<26>>>