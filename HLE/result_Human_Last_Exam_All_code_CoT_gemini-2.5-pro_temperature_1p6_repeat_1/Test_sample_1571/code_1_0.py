import sys
from io import StringIO

# Redirect stdout to capture prints
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

# --- Calculation Logic ---

# Initialize variables
total_life_loss = 0
equation_parts = []

print("Player A's turn proceeds as follows to maximize Player B's life loss:")

# Step 1: Pre-combat setup
print("\n1. Pre-Combat Phase:")
print("   - Player A taps one Mountain to cast Iron Apprentice for 1 mana.")
print("   - This provides a cheap artifact creature to sacrifice.")

# Step 2: Combat Phase
print("\n2. Combat Phase:")
# 2a. Attack triggers
clawing_torment_loss = 1
total_life_loss += clawing_torment_loss
equation_parts.append(str(clawing_torment_loss))
print(f"   - Player A attacks with all available creatures. The Clawing Torment attached to the attacking Scrap Welder causes Player B to lose {clawing_torment_loss} life.")

# 2b. Combat tricks and blocks
print("   - Player A activates Mukotai Soulripper's attack trigger, sacrificing the Iron Apprentice. Mukotai becomes a 6/6 creature with Menace.")
print("   - Player B must block the 6/6 Menace creature with both of their 2/1 Mothrider Patrols to minimize life loss. No damage from Mukotai gets through.")

# 2c. Combat damage
print("   - The remaining creatures are unblocked.")
ironhoof_boar_damage = 4
scrap_welder_damage = 3
replication_specialist_damage = 2

combat_damage = ironhoof_boar_damage + scrap_welder_damage + replication_specialist_damage
total_life_loss += combat_damage
equation_parts.append(str(ironhoof_boar_damage))
equation_parts.append(str(scrap_welder_damage))
equation_parts.append(str(replication_specialist_damage))
print(f"   - Ironhoof Boar deals {ironhoof_boar_damage} damage.")
print(f"   - Scrap Welder deals {scrap_welder_damage} damage.")
print(f"   - Replication Specialist deals {replication_specialist_damage} damage.")


# Step 3: Post-combat damage
print("\n3. Post-Combat Phase:")
# 3a. Twinshot Sniper Channel
twinshot_sniper_damage = 2
total_life_loss += twinshot_sniper_damage
equation_parts.append(str(twinshot_sniper_damage))
print(f"   - Player A uses their remaining mana to activate Twinshot Sniper's Channel ability from hand, dealing {twinshot_sniper_damage} damage to Player B.")

# 3b. March of Wretched Sorrow
march_damage = 2
total_life_loss += march_damage
equation_parts.append(str(march_damage))
print(f"   - Player A then casts March of Wretched Sorrow. They pay the {B} cost and exile the black card Junji from their hand, which adds {B}{B} to be spent on X. Thus X=2, and Player B loses {march_damage} life.")


# Step 4: Final Calculation
print("\n---")
print("Total life loss for Player B is the sum of all these sources:")
final_equation = " + ".join(equation_parts)
print(f"Final Calculation: {final_equation} = {total_life_loss}")

# --- Final Output ---
sys.stdout = old_stdout
output_str = captured_output.getvalue()
print(output_str)
print(f"<<<{total_life_loss}>>>")