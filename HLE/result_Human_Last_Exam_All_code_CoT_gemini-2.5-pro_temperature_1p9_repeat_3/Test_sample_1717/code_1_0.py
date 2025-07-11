# The base health of the Eye of Cthulhu in Terraria's Expert Mode for one player.
base_health_expert = 3640

# In Expert Mode, boss health increases for each additional player.
# This bonus health is capped at a maximum of 1000% (or a 10x bonus).
# The total maximum health multiplier is the base 100% + the 1000% bonus.
max_multiplier = 1.0 + 10.0

# Calculate the theoretical maximum health with an infinite number of players (which triggers the cap).
max_health = base_health_expert * max_multiplier

# The final equation to find the maximum health.
first_number = base_health_expert
second_number = max_multiplier

print(f"The theoretical maximum health of the Eye of Cthulhu is calculated as follows:")
print(f"{first_number} * {second_number} = {int(max_health)}")