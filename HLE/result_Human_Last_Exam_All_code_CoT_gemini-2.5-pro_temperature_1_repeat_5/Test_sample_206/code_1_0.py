# Based on the analysis of the board state, the goal is to deal 3 damage for a win.
# An all-out attack is the optimal strategy.

# 1. Define the attacking creatures.
# An attack with all creatures guarantees lethal damage.
attacking_creatures = [1, 2, 3, 4]

# 2. Determine which creatures will die from the optimal attack and blocks.
# - My (4) Wind Drake (2/2, Flying) attacks.
# - Opponent's (8) Skywinder Drake (2/2, Flying) must block it. They trade and both die.
# - My (1) Axegrinder Giant, (2) Centaur Courser, and (3) Axebane Beast attack.
# - Opponent's (7) Rusted Sentinel blocks the Giant and dies.
# - Opponent's (6) Grizzly Bears blocks the Courser and dies.
# - My (3) Axebane Beast is unblocked, dealing lethal damage.
dying_creatures = [4, 6, 7, 8]

# 3. Sort the lists in increasing order as requested.
attacking_creatures.sort()
dying_creatures.sort()

# 4. Format the output string.
# The format is "(c1), (c2); (d1), (d2)".
# We use a loop to build the string parts for attackers and dying creatures.
attackers_str = ", ".join(f"({c})" for c in attacking_creatures)
dying_str = ", ".join(f"({c})" for c in dying_creatures)

# 5. Print the final answer.
print(f"{attackers_str}; {dying_str}")