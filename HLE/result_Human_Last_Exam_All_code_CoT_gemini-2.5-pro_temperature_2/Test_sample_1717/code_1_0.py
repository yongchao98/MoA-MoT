# Step 1: Define the base health for the Eye of Cthulhu in Expert Mode for one player.
base_expert_health = 3640

# Step 2: Define the maximum health scaling multiplier for bosses in multiplayer.
# This value is a hard-coded limit in the game.
max_multiplier = 1000

# Step 3: Calculate the theoretical maximum health.
max_health = base_expert_health * max_multiplier

# Step 4: Print the final equation and the result.
# The format will show each number used in the calculation.
print(f"{base_expert_health} * {max_multiplier} = {max_health}")