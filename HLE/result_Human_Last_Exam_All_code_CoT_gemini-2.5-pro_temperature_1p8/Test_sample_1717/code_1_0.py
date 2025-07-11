# 1. Define the core variables based on the game's data.
base_health = 3640
player_multiplier = 0.35
max_players = 255

# 2. Calculate the health for the additional players.
# The scaling applies to all players after the first one.
additional_players = max_players - 1
health_scaling_factor = 1 + (player_multiplier * additional_players)

# 3. Calculate the final theoretical maximum health.
# The result is cast to an integer because health points in Terraria are whole numbers.
max_health = int(base_health * health_scaling_factor)

# 4. Print the final equation and the result.
print("The formula for max health is: Base Health * (1 + Multiplier * (Max Players - 1))")
print(f"The calculation is: {base_health} * (1 + {player_multiplier} * ({max_players} - 1)) = {max_health}")
print(f"The theoretical maximum health for the Eye of Cthulhu with {max_players} players in Expert Mode is: {max_health}")