# In Terraria, boss health scales with the number of players in Expert Mode.
# However, this scaling is not infinite and has a maximum cap.

# Step 1: Define the base health for the Eye of Cthulhu in single-player Expert Mode.
# This value is 3640 HP.
base_expert_health = 3640

# Step 2: Define the maximum health multiplier.
# The health scaling is capped at a 1000% increase, which means the final health
# cannot exceed 11 times the base single-player Expert Mode health.
# This cap is reached at 256 players, so for any number of players at or above 256
# (including a theoretical infinite number), this multiplier is used.
max_health_multiplier = 11

# Step 3: Calculate the theoretical maximum health.
max_theoretical_health = base_expert_health * max_health_multiplier

# Step 4: Print the explanation and the final equation.
print("The base health of the Eye of Cthulhu in Expert Mode is 3640.")
print("The health multiplier for players is capped at a maximum of 11x.")
print("The calculation for the theoretical maximum health is:")
print(f"{base_expert_health} * {max_health_multiplier} = {max_theoretical_health}")
