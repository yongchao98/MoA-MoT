import math

# Step 1: Explain the health scaling formula for the Eye of Cthulhu in Expert Mode.
# The formula is: Health = BaseHealth * (1 + ScalingFactor * (NumberOfPlayers - 1))
base_health_expert = 3640
scaling_factor_per_player = 0.35

print("The health of the Eye of Cthulhu in Expert Mode is calculated based on its base health and the number of players.")
print(f"Base Expert Health: {base_health_expert}")
print(f"Health increase per additional player: {scaling_factor_per_player * 100}% of base expert health")
print(f"Formula: Health = {base_health_expert} * (1 + {scaling_factor_per_player} * (NumberOfPlayers - 1))\n")

# Step 2: Explain the implication of infinite players.
print("As the 'NumberOfPlayers' approaches infinity, the calculated health also approaches infinity because the formula is linear.")
print("However, the game's source code stores NPC health in a 32-bit signed integer.\n")

# Step 3: Explain and calculate the technical limit.
print("This data type has a maximum possible value, which acts as the theoretical health cap for any NPC in the game.")
print("The calculation for this maximum value is 2 to the power of 31, minus 1.\n")

# Step 4: Print the final equation and the result.
max_health = 2**31 - 1
power = 31
base = 2

print("Final Equation and Theoretical Maximum Health:")
# The prompt requires printing each number in the final equation.
print(f"{base} ** {power} - 1 = {max_health}")