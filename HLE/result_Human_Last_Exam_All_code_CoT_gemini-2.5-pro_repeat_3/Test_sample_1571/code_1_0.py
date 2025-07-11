# Player A's unblocked creatures' power
junji_damage = 5
mukotai_soulripper_damage = 2
replication_specialist_damage = 2

# Player B's life loss is the sum of the damage from unblocked creatures.
# Player B's optimal blocking strategy is to block Ironhoof Boar (4 power)
# and Scrap Welder (3 power) to minimize the total damage taken.
# This leaves Junji, Mukotai Soulripper, and Replication Specialist unblocked.
total_life_lost = junji_damage + mukotai_soulripper_damage + replication_specialist_damage

print("To maximize damage, Player A casts Junji, the Midnight Sky and attacks with all available creatures.")
print("Player B's optimal defense is to block the Ironhoof Boar and Scrap Welder.")
print("The unblocked creatures are Junji, Mukotai Soulripper, and Replication Specialist.")
print(f"The total life Player B will lose is calculated by the sum of their power:")
print(f"{junji_damage} (from Junji) + {mukotai_soulripper_damage} (from Mukotai Soulripper) + {replication_specialist_damage} (from Replication Specialist) = {total_life_lost}")
