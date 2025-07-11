import collections

def calculate_tiling_ways(n):
    """
    Calculates T_n, the number of ways to tile a 2xn board using
    2x1, 2x2, and 2x4 tiles.
    """
    if n < 0:
        return 0
    
    # Memoization cache
    memo = collections.defaultdict(int)

    # T_0 = 1: One way to tile a 2x0 board (the empty tiling).
    memo[0] = 1

    # T_1 = 1: One vertical 2x1 tile.
    if n >= 1:
        memo[1] = 1

    # T_2 = 3: Two vertical 2x1s, two horizontal 2x1s, or one 2x2.
    if n >= 2:
        memo[2] = 3
        
    # T_3 = 5: (VVV, HHV, VHH) + ([2x2][V], [V][2x2])
    if n >= 3:
        memo[3] = 5

    # Use the recurrence T_n = T_{n-1} + 2*T_{n-2} + T_{n-4} for n >= 4
    # The recurrence can also generate T_2 and T_3 if T_0 and T_1 are set.
    # T_2 = T_1 + 2*T_0 = 1 + 2*1 = 3
    # T_3 = T_2 + 2*T_1 = 3 + 2*1 = 5
    # So we can run a loop from 2 to n.
    for i in range(2, n + 1):
        # We need T_{i-4}, which can be 0 if i-4 is negative.
        t_i_minus_4 = memo[i - 4] if i >= 4 else 0
        memo[i] = memo[i - 1] + 2 * memo[i - 2] + t_i_minus_4
        
    return memo

# We need to calculate T_4
target_n = 4
t_values = calculate_tiling_ways(target_n)

t0 = t_values[0]
t2 = t_values[2]
t3 = t_values[3]
t4 = t_values[4]

# The recurrence for T_4 is T_4 = T_3 + 2*T_2 + T_0
print(f"To calculate T_4, we use the recurrence relation T_n = T_(n-1) + 2*T_(n-2) + T_(n-4).")
print(f"For n=4, this is: T_4 = T_3 + 2*T_2 + T_0")
print()
print(f"First, we calculate the required base cases:")
print(f"T_0 = {t0}")
print(f"T_2 = {t2}")
print(f"T_3 = {t3}")
print()
print(f"Now, we substitute these values into the equation:")
print(f"T_4 = {t3} + 2 * {t2} + {t0} = {t4}")