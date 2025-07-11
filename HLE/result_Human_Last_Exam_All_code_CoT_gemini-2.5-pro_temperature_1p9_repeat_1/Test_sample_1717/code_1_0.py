import math

# Step 1: Establish the Health Scaling Formula
# According to the official Terraria Wiki, the Eye of Cthulhu's health in Expert Mode
# for a single player is 3640. For each additional player, health increases by 35%.
# The formula is: Health(N) = BaseExpertHealth * (1 + 0.35 * (N - 1))
# where N is the number of players.

# Let's simplify the formula:
# Health(N) = 3640 * (1 + 0.35*N - 0.35)
# Health(N) = 3640 * (0.65 + 0.35*N)
# Health(N) = 2366 + 1274*N

# Step 2: Analyze the Formula for Infinite Players
# The formula Health(N) = 2366 + 1274*N is a linear function.
# As the number of players N approaches infinity, the health value also approaches infinity.

# Step 3: Identify the Technical Limitation
# Given the mathematical formula tends to infinity, the theoretical maximum health is limited
# by the game's technical implementation. In Terraria's source code, an NPC's health
# (lifeMax) is stored as a 32-bit signed integer.

# Step 4: Calculate the Maximum Value
# The maximum value for a 32-bit signed integer is 2^31 - 1.

# Base for the exponentiation
base = 2
# Power for the exponentiation
power = 31
# Value to subtract
subtrahend = 1

# Final calculation
max_health = (base ** power) - subtrahend

print("The health of the Eye of Cthulhu in Expert Mode is a linear function of the number of players (N):")
print("Health = 2366 + 1274 * N")
print("\nAs N approaches infinity, the health is limited by the maximum value of a 32-bit signed integer.")
print("This technical limit is calculated as follows:")
# The format f-string with commas makes the large number readable.
print(f"({base}^{power}) - {subtrahend} = {max_health:,}")
