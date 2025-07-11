# Board dimensions and monster count
rows = 2024
columns = 2023
monsters = 2022

# Let C be the number of columns.
C = columns

# Let M be the number of monsters.
M = monsters

# The number of monsters is C - 1.
# This means there is exactly one "safe column" with no monsters.
# To guarantee a win, Turbo must identify this safe column.
# The worst-case scenario involves eliminating all unsafe columns until only one option remains.

# Turbo can probe the grid. Each failed attempt reveals one monster's location,
# thus eliminating one column from the set of potentially safe columns.

# After k attempts, in the worst case, k unsafe columns are identified.
# There are C - k columns remaining that could be the safe one.

# Let's consider the state after k = C - 2 attempts.
# Number of attempts made = C - 2 = 2023 - 2 = 2021.
# In the worst-case, 2021 monsters have been found, each in a different column.
# Number of columns whose status is unknown = C - (C - 2) = 2.
# One of these two columns is the safe one; the other contains the last monster.

# At the start of attempt k+1 = (C - 2) + 1 = C - 1, Turbo knows the
# final monster lies in one of two possible cells.
# He can devise a path that avoids BOTH of these potential monster locations
# by routing his path through one of the C-2 columns he already knows is unsafe,
# making a small detour around the specific cell where he knows the monster is.
# This guarantees success on this attempt.

# The number of attempts required is therefore (C - 2) probing attempts + 1 final successful attempt.
# Total attempts = (C - 2) + 1 = C - 1.

n = C - 1

print("The number of rows is R =", rows)
print("The number of columns is C =", columns)
print("The number of monsters is M =", monsters)
print("The number of safe columns is C - M =", columns - monsters)
print("The strategy requires determining the safe column.")
print("In the worst case, Turbo needs to eliminate C-2 unsafe columns to reduce the possibilities to two.")
print("This requires C-2 attempts.")
print("The (C-1)th attempt can then be a guaranteed success.")
print(f"The minimum value of n is C - 1 = {C} - 1 = {n}")