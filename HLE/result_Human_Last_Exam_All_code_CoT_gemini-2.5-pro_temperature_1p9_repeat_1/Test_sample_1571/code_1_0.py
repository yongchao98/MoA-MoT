import sys
# Redirect print to a string to control final output format.
from io import StringIO
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()

# --- Calculation of Optimal Play ---

# Initial state and variables
total_life_loss = 0
unblocked_damage_sources = {}

# Step 1: Channel Twinshot Sniper
channel_damage = 2
total_life_loss += channel_damage

# Step 2: Set up combat attackers
# After casting Junji and sacrificing a Ninja for Mukotai's ability
attackers = {
    "Mukotai Soulripper (buffed)": 4,
    "Ironhoof Boar": 4,
    "Scrap Welder": 3,
    "Replication Specialist": 2,
    "Ninja Token": 2
}

# Step 3: Player B blocks the 5/5 Junji, so it deals no damage to Player B.
# Calculate combat damage from unblocked attackers.
combat_damage = 0
for name, power in attackers.items():
    unblocked_damage_sources[name] = power
    combat_damage += power

total_life_loss += combat_damage


# --- Final Output Generation ---

# Reset stdout
sys.stdout = old_stdout

print("To calculate the maximum life loss for Player B, we follow Player A's optimal turn:")

print(f"\n1. First, Player A uses Twinshot Sniper's Channel ability. This is direct damage.")
print(f"   - Damage from Channel: {channel_damage}")


print(f"\n2. Next, Player A attacks. Player B's two 1/1 flying blockers must block the 5/5 flying, menace Junji, so Junji deals no damage to Player B.")
print("   The remaining creatures are unblocked and deal combat damage:")

equation_parts = [str(channel_damage)]
for source, damage in unblocked_damage_sources.items():
    print(f"   - Damage from {source}: {damage}")
    equation_parts.append(str(damage))

print("\n3. The total life loss is the sum of the Channel damage and the combat damage.")

final_equation = " + ".join(equation_parts)
print(f"The final calculation is:")
print(f"{final_equation} = {total_life_loss}")

<<<17>>>