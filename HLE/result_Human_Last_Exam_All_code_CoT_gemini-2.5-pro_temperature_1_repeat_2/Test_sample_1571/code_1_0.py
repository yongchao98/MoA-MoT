import collections

# --- State Initialization ---
# Player A's board state before playing cards from hand
attackers_power = {
    "Replication Specialist": 2,
    "Mukotai Soulripper": 2,
    "Scrap Welder": 3,
    "Ironhoof Boar": 4,
    "Ninja token 1": 1,
    "Ninja token 2": 1,
}
initial_artifact_count = 4

# --- Player A's Optimal Pre-Combat Plays ---
# Player A plays Iron Apprentice, replicates it, and plays Twinshot Sniper.
# This adds 3 artifacts to the board.
final_artifact_count = initial_artifact_count + 3
# It also adds two new creatures to the board for the attack.
attackers_power["Twinshot Sniper"] = 2
attackers_power["Iron Apprentice"] = 1 # The other is used to crew and is sacrificed

# --- Combat Phase ---
# 1. Resolve attack triggers
# Ironhoof Boar gets +X/+0 for X artifacts
attackers_power["Ironhoof Boar"] += final_artifact_count
# Mukotai Soulripper gets +2/+0 for sacrificing a creature
attackers_power["Mukotai Soulripper"] += 2

# 2. Determine Player B's optimal blocks
# Player B has two blockers and wants to minimize damage. They block the two
# highest-power creatures they can legally block. Mukotai has menace, so they
# block the two highest non-menace creatures.
# We create a list of non-menace attackers to find the highest power ones.
non_menace_attackers = {k: v for k, v in attackers_power.items() if k != "Mukotai Soulripper"}
sorted_non_menace_attackers = sorted(non_menace_attackers.items(), key=lambda item: item[1], reverse=True)

# Player B blocks the top two.
blocked_creature_1_name = sorted_non_menace_attackers[0][0]
blocked_creature_2_name = sorted_non_menace_attackers[1][0]

# 3. Identify unblocked creatures and calculate damage
unblocked_attackers = attackers_power.copy()
del unblocked_attackers[blocked_creature_1_name]
del unblocked_attackers[blocked_creature_2_name]

total_damage = sum(unblocked_attackers.values())

# --- Output the result ---
print(f"Player A's optimal play results in {final_artifact_count} artifacts.")
print("During combat, the attacker powers are:")
for name, power in sorted(attackers_power.items(), key=lambda item: item[1], reverse=True):
    print(f"- {name}: {power} power")

print(f"\nPlayer B minimizes damage by blocking '{blocked_creature_1_name}' ({attackers_power[blocked_creature_1_name]} power) and '{blocked_creature_2_name}' ({attackers_power[blocked_creature_2_name]} power).")

print("\nThe unblocked creatures deal damage equal to their power:")
# Sort for consistent output order
sorted_unblocked = collections.OrderedDict(sorted(unblocked_attackers.items(), key=lambda item: item[1], reverse=True))
equation_parts = [f"{power}" for power in sorted_unblocked.values()]
name_parts = [f"{name} ({power})" for name, power in sorted_unblocked.items()]
print(f"Unblocked: {', '.join(name_parts)}")
print(f"\nDamage Calculation: {' + '.join(equation_parts)} = {total_damage}")
print(f"The total amount of life Player B will lose is {total_damage}.")
