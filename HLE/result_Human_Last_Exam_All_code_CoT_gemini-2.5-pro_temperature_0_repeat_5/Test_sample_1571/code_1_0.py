# Define the power of Player A's creatures involved in the final attack.
junji_power = 5
soulripper_base_power = 2
boar_power = 4
welder_power = 3
ninja_power = 1
apprentice_token_power = 1

# Calculate the final power of Mukotai Soulripper after its abilities resolve.
# It gets +2/+0 from its own ability and a +1/+1 counter from the sacrificed Apprentice.
soulripper_final_power = soulripper_base_power + 2 + 1

# Player B's optimal strategy is to use both flying blockers to block the flying menace Junji.
# This prevents damage from Junji, but all other creatures are unblocked.
# The total damage is the sum of the power of all unblocked creatures.
total_damage = soulripper_final_power + boar_power + welder_power + ninja_power + ninja_power + apprentice_token_power

# Print the final calculation step-by-step.
print("Player B minimizes damage by using both flying blockers on the 5/5 flying menace Junji.")
print("The remaining unblocked creatures will deal damage equal to their total power.")
print("The final damage calculation is:")
print(f"{soulripper_final_power} (from Mukotai Soulripper) + {boar_power} (from Ironhoof Boar) + {welder_power} (from Scrap Welder) + {ninja_power} (from Ninja) + {ninja_power} (from Ninja) + {apprentice_token_power} (from Apprentice Token) = {total_damage}")
