def calculate_tilings(n_max):
    """
    Calculates T_n for n from 0 to n_max using the recurrence relation.
    T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
    """
    T = {}
    # Base case for the recurrence calculation, T_n = 0 for n < 0
    # T_0 = 1 is the actual first term.
    
    for i in range(n_max + 1):
        if i == 0:
            T[i] = 1
            continue
        
        t_i_minus_1 = T.get(i - 1, 0)
        t_i_minus_2 = T.get(i - 2, 0)
        t_i_minus_4 = T.get(i - 4, 0)
        
        T[i] = t_i_minus_1 + 2 * t_i_minus_2 + t_i_minus_4
        
    return T

# Calculate tilings up to T_4
n = 4
T = calculate_tilings(n)

# Get the specific values needed for the T_4 calculation
t4 = T[4]
t3 = T[3]
t2 = T[2]
t0 = T[0]

# Print the explanation and the final equation
print(f"To find T_4, we use the recurrence relation T_n = T_(n-1) + 2*T_(n-2) + T_(n-4).")
print(f"For n=4, the equation is T_4 = T_3 + 2 * T_2 + T_0.")
print(f"We have the following values:")
print(f"T_3 = {t3}")
print(f"T_2 = {t2}")
print(f"T_0 = {t0}")
print(f"\nPlugging these values into the equation:")
print(f"T_4 = {t3} + 2 * {t2} + {t0} = {t4}")