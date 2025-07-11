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

# VP calculations
# 1. Alliance VPs
chaeilki_vp = developed_planets * 1  # 1 VP per developed planet
humans_vp = colonized_planets * 1   # 1 VP per colonized planet
usud_vp = conquered_planets * 1     # 1 VP per conquered planet

# 2. Ideology VPs
legarchaea_vp = prestige // 3  # 1 VP per 3 prestige

# 3. General VPs
credits_vp = credits // 5
productivity_vp = productivity // 5
discovery_vp = discovery // 5
influence_vp = influence // 5
relics_vp = relics * 1
cruisers_vp = cruisers // 2

# Total VP
total_vp = (chaeilki_vp + humans_vp + usud_vp + legarchaea_vp +
            credits_vp + productivity_vp + discovery_vp + influence_vp +
            relics_vp + cruisers_vp)

# Print the detailed scoring breakdown and final result
print("Calculating Victory Points...\n")
print(f"Alliance VPs:")
print(f"  Chaeilki (1 per developed planet): {chaeilki_vp}")
print(f"  Humans (1 per colonized planet): {humans_vp}")
print(f"  Us'ud (1 per conquered planet): {usud_vp}")
print("\nIdeology VPs:")
print(f"  Legarchaea (1 per 3 prestige): {legarchaea_vp}")
print("\nGeneral VPs:")
print(f"  Credits (1 per 5): {credits_vp}")
print(f"  Productivity (1 per 5): {productivity_vp}")
print(f"  Discovery (1 per 5): {discovery_vp}")
print(f"  Influence (1 per 5): {influence_vp}")
print(f"  Relics (1 per relic): {relics_vp}")
print(f"  Cruisers (1 per 2): {cruisers_vp}")
print("\n--------------------")
print("Final VP Equation:")
print(f"{chaeilki_vp} (Chaeilki) + {humans_vp} (Humans) + {usud_vp} (Us'ud) + "
      f"{legarchaea_vp} (Legarchaea) + {credits_vp} (Credits) + {productivity_vp} (Productivity) + "
      f"{discovery_vp} (Discovery) + {influence_vp} (Influence) + {relics_vp} (Relics) + {cruisers_vp} (Cruisers)")
print("--------------------")
print(f"Total VP: {total_vp}")
