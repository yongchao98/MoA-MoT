# Player A aims to maximize the life Player B loses.
# Player B aims to minimize the life they lose.

# Initial state:
# Player B has two 1/1 flying blockers and no ground blockers.
# Player A has several ground creatures and a hand with direct damage potential.

life_lost_components = []
life_lost_reasons = []

# 1. Pre-combat Main Phase: Player A uses Twinshot Sniper's Channel ability.
# This is the most efficient use of mana for direct damage.
channel_damage = 2
life_lost_components.append(channel_damage)
life_lost_reasons.append(f"Player A uses Twinshot Sniper's Channel ability, dealing {channel_damage} damage.")

# 2. Combat Phase: Declare Attackers and Resolve Triggers.
# Player A crews Mukotai Soulripper by tapping the two 1/1 Ninja tokens.
# Player A attacks with all available ground creatures.
# The Clawing Torment on Scrap Welder triggers upon attacking.
trigger_life_loss = 1
life_lost_components.append(trigger_life_loss)
life_lost_reasons.append(f"Scrap Welder attacks, and the ability from Clawing Torment causes Player B to lose {trigger_life_loss} life.")

# 3. Combat Damage: Player B cannot block any attackers.
# Damage from Mukotai Soulripper (4/5)
mukotai_damage = 4
life_lost_components.append(mukotai_damage)
life_lost_reasons.append(f"Mukotai Soulripper deals {mukotai_damage} unblocked combat damage.")

# Damage from Scrap Welder (2/2)
scrap_welder_damage = 2
life_lost_components.append(scrap_welder_damage)
life_lost_reasons.append(f"Scrap Welder deals {scrap_welder_damage} unblocked combat damage.")

# Damage from Ironhoof Boar (4/3)
ironhoof_boar_damage = 4
life_lost_components.append(ironhoof_boar_damage)
life_lost_reasons.append(f"Ironhoof Boar deals {ironhoof_boar_damage} unblocked combat damage.")

# Damage from Replication Specialist (3/4)
replication_specialist_damage = 3
life_lost_components.append(replication_specialist_damage)
life_lost_reasons.append(f"Replication Specialist deals {replication_specialist_damage} unblocked combat damage.")

# 4. Final Calculation
total_life_lost = sum(life_lost_components)

# Print the breakdown of life loss
print("Here is the breakdown of how much life Player B loses:")
for reason in life_lost_reasons:
    print(f"- {reason}")

# Print the final equation and total
equation_str = " + ".join(map(str, life_lost_components))
print("\nThe total life lost by Player B is calculated as follows:")
print(f"Total Life Lost = {equation_str} = {total_life_lost}")
