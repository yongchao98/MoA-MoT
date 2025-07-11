# Define the constants from the Terraria health scaling formula for Expert Mode.
base_health = 3640
health_increase_per_player = 0.35

# The formula for total health is:
# Total Health = base_health * (1 + (number_of_players - 1) * health_increase_per_player)

print("The formula for the Eye of Cthulhu's health in Expert Mode is based on the number of players.")
print("Final Health = Base Health * (1 + (Number of Players - 1) * Health Increase per Player)")
print("\nHere are the numbers in that final equation:")
print(f"Base Health: {base_health}")
print(f"Constant '1': 1")
print(f"Constant '-1': -1")
print(f"Health Increase per Player: {health_increase_per_player}")

print("\nSince the 'Number of Players' can theoretically be infinite, and the health scales linearly with it,")
print("the theoretical maximum health is also infinite. There is no cap on this scaling in the game's code.")
