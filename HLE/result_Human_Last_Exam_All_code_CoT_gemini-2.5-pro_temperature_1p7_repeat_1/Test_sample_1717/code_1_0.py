# Define the base values for the calculation based on Terraria's Expert Mode mechanics.
base_health_expert = 3640  # Health for the Eye of Cthulhu with one player.
scaling_factor_expert = 0.35 # Health increases by 35% for each additional player.

# Calculate the amount of health added for each player after the first.
health_increase_per_player = int(base_health_expert * scaling_factor_expert)

print("To determine the theoretical maximum health of the Eye of Cthulhu, we use the multiplayer scaling formula from Terraria's source code.")
print("-" * 50)
print(f"1. Base Health in Expert Mode (1 player): {base_health_expert}")
print(f"2. Health Increase per Additional Player: {scaling_factor_expert * 100}% of base health")
print(f"   Calculation: {base_health_expert} * {scaling_factor_expert} = {health_increase_per_player} health per extra player.")
print("-" * 50)
print("The general formula for health (H) based on the number of players (n) is:")
print("H(n) = BaseHealth + (n - 1) * HealthIncreasePerPlayer")
print("\nPlugging in the values for the Eye of Cthulhu:")
print(f"Total Health = {base_health_expert} + (NumberOfPlayers - 1) * {health_increase_per_player}")
print("-" * 50)
print("\nConclusion for an infinite number of players:")
print("The health formula is a linear equation. As the 'NumberOfPlayers' increases, the 'Total Health' increases without any in-game limit.")
print("Therefore, with a theoretical infinite number of players, the Eye of Cthulhu's health would also approach infinity.")
