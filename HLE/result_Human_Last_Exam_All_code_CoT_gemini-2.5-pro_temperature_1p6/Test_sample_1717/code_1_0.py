import sys

# According to the official Terraria Wiki, which references the game's source code,
# the health of bosses in Expert Mode scales with the number of players.

# Step 1: Define the base health for the Eye of Cthulhu in Expert Mode.
base_health_expert = 3640

# Step 2: Define the health scaling factor for each player after the first.
# Health increases by 35% for each additional player.
per_player_increase = 0.35

# Step 3: Present the full formula.
# Let 'N' be the number of players. The multiplier for health is (1 + 0.35 * (N - 1)).
# The '1' represents the base health for the first player.
# The '0.35 * (N - 1)' represents the additional 35% health for every player beyond the first.
print("The formula for the Eye of Cthulhu's health in Expert Mode is based on the number of players (N):")
print("\nTotal Health = Base Health * (1 + Per-Player Increase * (N - 1))")
print("\nPlugging in the numbers, the equation is:")
print(f"Total Health = {base_health_expert} * (1 + {per_player_increase} * (N - 1))")
print("-" * 50)

# Step 4: Analyze the formula for an infinite number of players.
print("The question asks for the theoretical maximum health given an infinite number of players.")
print("\nThis means we need to evaluate the formula as 'N' approaches infinity.")
print("\n1. As the number of players 'N' approaches infinity, the term '(N - 1)' also approaches infinity.")
print(f"\n2. As '(N - 1)' approaches infinity, the health multiplier '(1 + {per_player_increase} * (N - 1))' also approaches infinity.")
print(f"\n3. When a positive constant ({base_health_expert}) is multiplied by a value that approaches infinity, the result also approaches infinity.")
print("-" * 50)

# Final Conclusion
print("\nConclusion:")
print("The health of the Eye of Cthulhu scales linearly with the number of players and has no upper limit or cap.")
print("Therefore, with a theoretical infinite number of players, the Eye of Cthulhu would have a theoretical infinite amount of health.")
