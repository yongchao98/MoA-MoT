# Define the constants based on Terraria's game mechanics for Expert Mode.
base_health = 3640
health_increase_per_player = 0.35

# The problem asks for the health with an infinite number of players.
# The health formula is a linear function: H(n) = base_health * (1 + 0.35 * (n - 1)),
# where 'n' is the number of players. As 'n' approaches infinity, H(n) also approaches infinity.
# Therefore, there is no finite maximum health.

print("The theoretical maximum health of the Eye of Cthulhu in Expert Mode with an infinite number of players is infinite.")
print("This is because its health scales linearly with the number of players.")
print("Formula: Total Health = BaseHealth * (1 + 0.35 * (NumberOfPlayers - 1))\n")

# To demonstrate the formula, we will calculate the health for a very large,
# but finite, number of players, for example: 1,000,000.
number_of_players = 1000000

print(f"--- Example Calculation for {number_of_players:,} Players ---")

# Calculate the final health using the formula.
# We calculate the scaling factor first for clarity.
scaling_factor = 1 + health_increase_per_player * (number_of_players - 1)
final_health = base_health * scaling_factor

# Print the final equation showing how each number is used, as requested.
# Health is an integer value in the game.
print("\nFinal Equation:")
print(f"{int(final_health)} = {base_health} * (1 + {health_increase_per_player} * ({number_of_players} - 1))")