# Step 1: Calculate the maximum damage from March of Wretched Sorrow.
# Player A has 6 mana and can exile one black card for a cost reduction of 2.
# The spell's cost is X + 1 (for the black mana pip).
# The total "payment" available is 6 mana + 2 from the cost reduction = 8.
# X + 1 = 8, so X = 7.
damage_from_march = 7

# Step 2: Define the power of unblocked creatures in combat.
# Player B must block the highest-power non-menace creatures.
# Attackers' powers: Ironhoof Boar (4), Replication Specialist (2),
# Mukotai Soulripper (2, menace), Scrap Welder (1), Ninja (1).
# B blocks the 4-power Boar and 2-power Specialist.
# Unblocked attackers are Mukotai, Scrap Welder, and the Ninja.
damage_from_mukotai = 2
damage_from_scrap_welder = 1
damage_from_ninja = 1

# Step 3: Calculate the total combat damage.
total_combat_damage = damage_from_mukotai + damage_from_scrap_welder + damage_from_ninja

# Step 4: Calculate the total life lost by Player B.
total_life_loss = damage_from_march + total_combat_damage

# Step 5: Print the final equation and the result.
print("Player B's total life loss is calculated by summing the damage from the spell and the unblocked creatures in combat.")
print(f"({damage_from_march}) + ({damage_from_mukotai}) + ({damage_from_scrap_welder}) + ({damage_from_ninja}) = {total_life_loss}")
